"""
zmp.py

Functions associated with zmp inversion
"""

from functools import reduce
import matplotlib.pyplot as plt

import math

import psi4
psi4.core.be_quiet()
import numpy as np
eps = np.finfo(float).eps
from opt_einsum import contract
from scipy.optimize import minimize

from psi4.driver.p4util.solvers import DIIS

import sys

class ZMP():

############### ZMP SCF

    def zmp_with_scf(self, 
            lam,
            opt_max_iter, 
            opt_tol, 
            ):

        """
        Performs ZMP optimization according to: 
        
        1) 'From electron densities to Kohn-Sham kinetic energies, orbital energies,
        exchange-correlation potentials, and exchange-correlation energies' by
        Zhao + Morrison + Parr. 
        https://doi.org/10.1103/PhysRevA.50.2138

        Additional DIIS algorithms obtained from:
        2) 'Psi4NumPy: An interactive quantum chemistry programming environment 
        for reference implementations and rapid development.' by 
        Daniel G.A. Smith and others. 
        https://doi.org/10.1021/acs.jctc.8b00286
        """

        self.shift = 0.1 * lam
        self.diis_space = 100

        print("Running SCF ZMP:")
        self.zmp_scf(lam, opt_max_iter, D_conv=opt_tol)

    def zmp_scf(self, 
            lam = 100,
            maxiter = 100, 
            D_conv = psi4.core.get_option("SCF", "D_CONVERGENCE")):
        """
        Performs scf cycle
        Parameters
        ----------
        lam: integer, opt
            Global lagrange multiplier for effective potential that drives SCF calculation. 
            Defined in equation 7 and 8 of reference (1).
        maxiter: integer, opt
            Maximum number of iterations for SCF calculation
        D_conv: float
            Convergence parameter for density and DIIS error. 
            Default is Psi4's Density convergence parameter: 1e-06
        """

        A = np.array(self.A)
        S2 = self.S2
        H = self.T + self.V 

        # Trial & Residual Vector Lists
        state_vectors_a, state_vectors_b = [], []
        error_vectors_a, error_vectors_b = [], []

        #"Initial Guess for SCF
        Cocca = self.ct[0]
        Coccb = self.ct[1]
        Da = self.nt[0]
        Db = self.nt[1]
        D_old = Da.copy()

        # Obtain slice of target density to test convergence on grid
        # at the end of scf. 
        x = np.linspace(-10,10,201)
        y = np.linspace(-10,10,201)
        z = np.linspace(-10,10,201)
        grid = np.concatenate((x[:,None], y[:,None], z[:,None]), axis=1).T
        density0 = self.on_grid_density(grid, Da=self.nt[0], Db=self.nt[1])

        dd_old = 0.0

        print("-------------SCF TIME ----------------------")
        for SCF_ITER in range(1,maxiter):

#------------->  Generate Fock Matrix:
            Fa = np.zeros((self.nbf, self.nbf))
            Fb = np.zeros((self.nbf, self.nbf))
            Ja = contract('pqrs,rs->pq', self.I, Da)
            Jb = contract('pqrs,rs->pq', self.I, Db)

            #Equation 7 of Reference (1)
            v_c = lam * ( (Ja + Jb) - (self.J0[0] + self.J0[1]) )

            #Equation 10 of Reference (1)
            Fa += H + self.va + v_c 
            Fb += H + self.vb + v_c

            #Level Shift
            Fa += (S2 - reduce(np.dot, (S2, Da, S2)) * self.shift)
            Fb += (S2 - reduce(np.dot, (S2, Db, S2)) * self.shift)

#------------->  DIIS:
            if SCF_ITER > 1:
                #Construct the AO gradient
                # r = (A(FDS - SDF)A)_mu_nu
                grad_a = A.dot(Fa.dot(Da).dot(S2) - S2.dot(Da).dot(Fa)).dot(A)
                grad_a[np.abs(grad_a) < eps] = 0.0

                if SCF_ITER -1 < self.diis_space:
                    state_vectors_a.append(Fa.copy())
                    error_vectors_a.append(grad_a.copy())
                else:
                    state_vectors_a.append(Fa.copy())
                    error_vectors_a.append(grad_a.copy())

                #Build inner product of error vectors
                Bdim = len(state_vectors_a)
                Ba = np.empty((Bdim + 1, Bdim + 1))
                Ba[-1, :] = -1
                Ba[:, -1] = -1
                Ba[-1, -1] = 0
                Bb = Ba.copy()

                for i in range(len(state_vectors_a)):
                    for j in range(len(state_vectors_a)):
                        Ba[i,j] = np.einsum('ij,ij->', error_vectors_a[i], error_vectors_a[j])

                #Build almost zeros matrix to generate inverse. 
                RHS = np.zeros(Bdim + 1)
                RHS[-1] = -1

                #Find coefficient matrix:
                Ca = np.linalg.solve(Ba, RHS.copy())
                Ca[np.abs(Ca) < eps] = 0.0 

                #Generate new fock Matrix:
                Fa = np.zeros_like(Fa)
                for i in range(Ca.shape[0] - 1):
                    Fa += Ca[i] * state_vectors_a[i]

                if self.ref ==  1:
                    Fb = Fa.copy()
                    diis_error_a = np.max(np.abs(error_vectors_a[-1]))
                    diis_error = 2 * diis_error_a

                else:
                    grad_b = A.dot(Fa.dot(Db).dot(S2) - S2.dot(Da).dot(Fa)).dot(A)
                    grad_b[np.abs(grad_b) < eps] = 0.0
                    
                    if SCF_ITER -1 < self.diis_space:
                        state_vectors_b.append(Fb.copy())
                        error_vectors_b.append(grad_b.copy())
                    else:
                        state_vectors_b.append(Fb.copy())
                        error_vectors_b.append(grad_b.copy())

                    for i in range(len(state_vectors_a)):
                        for j in range(len(state_vectors_a)):
                            Bb[i,j] = np.einsum('ij,ij->', error_vectors_b[i], error_vectors_b[j])

                    diis_error_b = np.max(np.abs(error_vectors_b[-1]))
                    diis_error = diis_error_a + diis_error_b

                    Cb = np.linalg.solve(Bb, RHS.copy())
                    Cb[np.abs(Cb) < eps] = 0.0 

                    Fb = np.zeros_like(Fb)
                    for i in range(Cb.shape[0] - 1):
                        Fb += Cb[i] * state_vectors_b[i]

            else:
                diis_error = 1.0
            
#------------->  Diagonalization | Check convergence:

            Ca, Cocca, Da, eigs_a = self.diagonalize(Fa, self.nalpha)
            if self.ref == 2:
                Cb, Coccb, Db, eigs_b = self.diagonalize(Fb, self.nbeta)
            else: 
                Cb, Coccb, Db, eigs_b = Ca.copy(), Cocca.copy(), Da.copy(), eigs_a.copy()

            ddm = D_old - Da
            D_old = Da
            derror = np.max(np.abs(ddm))

            if True and np.mod(SCF_ITER,5) == 0.0:
                print(f"Iteration: {SCF_ITER:d} | Self Convergence Error: {derror:10.5e} | DIIS Error: {diis_error:10.5e}")
            if abs(derror) < D_conv and abs(diis_error) < D_conv:
                print("Done")
                break
            if SCF_ITER == maxiter - 1:
                raise ValueError("Maximum Number of SCF cycles reached. Try different settings.")

        density_current = self.on_grid_density(grid, Da=Da, Db=Db)
        grid_diff = np.max(np.abs(density0 - density_current))
        print(f"SCF Procedure successfull. Max density difference: {grid_diff}")
        self.Da = Da
        self.Db = Db

        #VXC is hartree-like Potential. We remove Fermi_Amaldi Guess. 
        self.proto_density_a = self.lam * (self.Da) - (self.lam + 1/(self.nalpha + self.nbeta)) * (self.nt[0])
        self.proto_density_b = self.lam * (self.Db) - (self.lam + 1/(self.nbeta + self.nalpha)) * (self.nt[1])


######## OPTIMIZATION ZMP

    def zmp(self, lam, opt_max_iter, opt_tol):
        """
        Calls scipy minimizer to minimize lagrangian. 
        """

        if self.ref == 2:
            raise ValueError("Unrestricted calculation has not been implemented. Try a Restricted calculation")

        self.v0 = np.zeros_like( self.T.flatten() )
        self.shift = 0.1 * lam

        self.Da = self.nt[0]
        self.Db = self.nt[1]

        opt_results = minimize( fun = self.lagrangian_zmp,
                                x0  = self.v0, 
                                jac = self.gradient_zmp,
                                # args=(diisa_obj, diisb_obj), 
                                method  = 'L-BFGS-B',
                                tol     = opt_tol,
                                options = {"maxiter" : opt_max_iter,
                                            "disp"    : False,}
                                )

        if opt_results.success == False:
            raise ValueError("Optimization was unsucessful, try a different intitial guess")
        else:

            npoints=201
            x = np.linspace(-6,6,npoints)[:,None]
            y = np.linspace(-6,6,npoints)[:,None]
            z = np.linspace(-6,6,npoints)[:,None] 
            grid = np.concatenate((x,y,z), axis=1).T
            density = self.on_grid_density(grid=grid, Da=self.Da, Db=self.Db)
            density0 = self.on_grid_density(grid=grid, Da=self.wfn.Da().np, Db=self.wfn.Db().np)

            dd = density[:,0] - density0[:,0]

            print("DD max innacuracy", np.max(dd))
            print("Optimization Successful")
            self.v_opt = opt_results.x
            self.opt_info = opt_results
            self.proto_density_a = self.lam * (self.Da) - (self.lam + 1/(self.nalpha + self.nbeta)) * (self.nt[0])
            self.proto_density_b = self.lam * (self.Db) - (self.lam + 1/(self.nalpha + self.nbeta)) * (self.nt[1])

    def diagonalize_with_guess_zmp(self, v): 
        """
        Diagonalize Fock matrix with additional external potential
        """
        v_c = self.lam * (2 * v) - self.lam * (np.array(self.J0[0]) + np.array(self.J0[1]))
        Fa = self.V + self.T + self.va + v_c
        Fa += (self.S2 - reduce(np.dot, (self.S2, self.Da, self.S2)) * self.shift)
        self.Ca, self.Coca, self.Da, self.eigvecs_a = self.diagonalize( Fa, self.nalpha )

        if self.ref == 1:
            self.Cb, self.Coca, self.Db, self.eigvecs_b = self.Ca.copy(), self.Coca.copy(), self.Da.copy(), self.eigvecs_a.copy()
        else:
            Fb = self.V + self.T + self.vb + v_c
            Fb += (self.S2 - reduce(np.dot, (self.S2, self.Db, self.S2)) * self.shift)
            self.Cb, self.Cocb, self.Db, self.eigvecs_b = self.diagonalize( Fb, self.nbeta )

    def lagrangian_zmp(self, v):
        """
        Lagrangian to be minimized wrt external potential
        Equation (5) of main reference
        """
        v = np.reshape(v, (self.nbf, self.nbf))
        self.diagonalize_with_guess_zmp(v)

        self.grad_a = self.Da - self.nt[0]
        self.grad_b = self.Db - self.nt[1] 

        kinetic    = np.sum(self.T * self.Da)
        potential  = np.sum((self.V + self.va) * (self.Da - self.nt[0]))
        optimizing = np.sum( v *  self.grad_a)

        if self.ref == 1:
            L = 2 * (kinetic + potential + optimizing)

        else:        
            kinetic    += np.sum((self.T * self.Db))
            potential  += np.sum((self.V + self.vb) * (self.Db - self.nt[1]))
            optimizing += np.sum( (v *  self.grad_b))

            L = kinetic + optimizing + potential
            # L = kinetic + optimizing

        if True:
            print(f"Kinetic: {kinetic:6.4f} | Potential: {np.abs(potential):6.4e} | From Optimization: {np.abs(optimizing):6.4e}")

        return -L

    def gradient_zmp(self, v):
        """
        Calculates gradient wrt target density
        Equation (11) of main reference
        """
        
        v = np.reshape(v, (self.nbf, self.nbf))
        self.diagonalize_with_guess_zmp(v)
        self.grad_a = self.Da - self.nt[0]
        self.grad_b = self.Db - self.nt[1] 

        if self.ref == 1:
            self.grad   = self.grad_a
        else:
            self.grad   = np.concatenate(( self.grad_a, self.grad_b ))

        self.grad = self.grad.flatten()
        return -self.grad
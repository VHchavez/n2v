"""
wuyang.py

Functions associated with wuyang inversion
"""

import numpy as np
from opt_einsum import contract
from scipy.optimize import minimize

class WuYang():
    """
    Performs Optimization as in: 10.1063/1.1535422 - Qin Wu + Weitao Yang

    Attributes:
    -----------
    lambda_rgl: {None, float}. If float, lambda-regularization is added with lambda=lambda_rgl.
    """

    regul_norm = None  # Regularization norm: ||v||^2
    lambda_reg = None  # Regularization constant

    def wuyang(self, opt_method, opt_max_iter, opt_tol):
        """
        Calls scipy minimizer to minimize lagrangian. 
        """

        if opt_method.lower() == 'bfgs' or opt_method.lower() == 'l-bfgs-b':
            opt_results = minimize( fun = self.lagrangian,
                                    x0  = self.v0, 
                                    jac = self.gradient,
                                    method  = opt_method,
                                    tol     = opt_tol,
                                    options = {"maxiter" : opt_max_iter,
                                                "disp"    : False,}
                                    )

        else:
            opt_results = minimize( fun = self.lagrangian,
                                    x0  = self.v0, 
                                    jac = self.gradient,
                                    hess = self.hessian,
                                    method = opt_method,
                                    tol    = opt_tol,
                                    options = {"maxiter"  : opt_max_iter,
                                                "disp"    : False, }
                                    )

        if opt_results.success == False:
            raise ValueError("Optimization was unsucessful (|grad|=%.2e) within %i iterations, "
                             "try a different intitial guess."% (np.linalg.norm(opt_results.jac), opt_results.nit)
                             + opt_results.message)
        else:
            print("Optimization Successful within %i iterations! "
                  "|grad|=%.2e" % (opt_results.nit, np.linalg.norm(opt_results.jac)))
            self.v_opt = opt_results.x
            self.opt_info = opt_results

        self.finalize_energy()

        # if debug=True:
        #     self.density_accuracy()

    def diagonalize_with_guess(self, v):
        """
        Diagonalize Fock matrix with additional external potential
        """
        vks_a = contract("ijk,k->ij", self.S3, v[:self.npbs]) + self.va
        fock_a = self.V + self.T + vks_a 
        self.Ca, self.Coca, self.Da, self.eigvecs_a = self.diagonalize( fock_a, self.nalpha )

        if self.ref == 1:
            self.Cb, self.Coca, self.Db, self.eigvecs_b = self.Ca.copy(), self.Coca.copy(), self.Da.copy(), self.eigvecs_a.copy()
        else:
            vks_b = contract("ijk,k->ij", self.S3, v[self.npbs:]) + self.vb
            fock_b = self.V + self.T + vks_b        
            self.Cb, self.Cocb, self.Db, self.eigvecs_b = self.diagonalize( fock_b, self.nbeta )

    def lagrangian(self, v):
        """
        Lagrangian to be minimized wrt external potential
        Equation (5) of main reference
        """

        self.diagonalize_with_guess(v)
        self.grad_a = contract('ij,ijt->t', (self.Da - self.nt[0]), self.S3)  
        self.grad_b = contract('ij,ijt->t', (self.Db - self.nt[1]), self.S3)

        kinetic     =   np.sum(self.T * (self.Da))
        potential   =   np.sum((self.V + self.va) * (self.Da - self.nt[0]))
        optimizing  =   np.sum(v[:self.npbs] * self.grad_a)

        if self.ref == 1:
            L = 2 * (kinetic + potential + optimizing)

        else:
            kinetic    +=   np.sum(self.T * (self.Db))
            potential  +=   np.sum((self.V + self.vb) * (self.Db - self.nt[1]))
            optimizing +=   np.sum(v[self.npbs:] * self.grad_b)
            L = kinetic + potential + optimizing

        # Add lambda-regularization
        if self.lambda_reg is not None:
            T = self.T_pbs
            if self.ref == 1:
                norm = 2 * (v[:self.npbs] @ T @ v[:self.npbs])
            else:
                norm = (v[self.npbs:] @ T @ v[self.npbs:]) + (v[:self.npbs] @ T @ v[:self.npbs])

            L -= norm * self.lambda_reg
            self.regul_norm = norm

        # if print_flag:
        #    print(f"Kinetic: {kinetic:6.4f} | Potential: {np.abs(potential):6.4e} | From Optimization: {np.abs(optimizing):6.4e}")

        return - L

    def gradient(self, v):
        """
        Calculates gradient wrt target density
        Equation (11) of main reference
        """
        self.diagonalize_with_guess(v)
        self.grad_a = contract('ij,ijt->t', (self.Da - self.nt[0]), self.S3)
        self.grad_b = contract('ij,ijt->t', (self.Db - self.nt[1]), self.S3) 

        if self.ref == 1:
            self.grad   = self.grad_a
        else:
            self.grad   = np.concatenate(( self.grad_a, self.grad_b ))

        if self.lambda_reg is not None:
            T = self.T_pbs
            if self.ref == 1:
                rgl_vector = 4 * self.lambda_reg*np.dot(T, v[:self.npbs])
                self.grad -= rgl_vector
            else:
                self.grad[:self.npbs] -= 2 * self.lambda_reg*np.dot(T, v[:self.npbs])
                self.grad[self.npbs:] -= 2 * self.lambda_reg*np.dot(T, v[self.npbs:])

        return -self.grad

    def hessian(self, v):
        """
        Calculates gradient wrt target density
        Equation (13) of main reference
        """

        self.diagonalize_with_guess(v)

        na, nb = self.nalpha, self.nbeta

        eigs_diff_a = self.eigvecs_a[:na, None] - self.eigvecs_a[None, na:]
        C3a = contract('mi,va,mvt->iat', self.Ca[:,:na], self.Ca[:,na:], self.S3)
        Ha = 2 * contract('iau,iat,ia->ut', C3a, C3a, eigs_diff_a**-1)

        if self. ref == 1:
            if self.lambda_reg is not None:
                Ha -= 4 * self.T_pbs * self.lambda_reg
            Hs = Ha

        else:

            eigs_diff_b = self.eigvecs_b[:nb, None] - self.eigvecs_b[None, nb:]
            C3b = contract('mi,va,mvt->iat', self.Cb[:,:nb], self.Cb[:,nb:], self.S3)
            Hb = 2 * contract('iau,iat,ia->ut', C3b, C3b, eigs_diff_b**-1)
            if self.lambda_reg is not None:
                Ha -= 2 * self.T_pbs * self.lambda_reg
                Hb -= 2 * self.T_pbs * self.lambda_reg
            Hs = np.block(
                            [[Ha,                               np.zeros((self.npbs, self.npbs))],
                            [np.zeros((self.npbs, self.npbs)), Hb                              ]]
                        )

        return - Hs

    def find_regularization_constant(self, guide_potential_components, opt_method="trust-krylov",
                                     tol=None, opt=None, lambda_list=None):
        """
        Finding regularization constant lambda.

        Note: it is recommend to set a specific convergence criteria by opt or tol,
                in order to control the same convergence
                for different lambda value.

        After the calculation is done, one can plot the returns to select a good lambda.

        Parameters:
        -----------
        guide_potential_components: a list of string
            the components for guide potential v0.
            see Inverter.generate_components() for details.

        opt_method: string default: "trust-krylov"
            opt_methods available in scipy.optimize.minimize

        tol: float
            Tolerance for termination. See scipy.optimize.minimize for details.

        opt: dictionary, optional
            if given:
                scipy.optimize.minimize(method=opt_method, options=opt)

        lambda_list: np.ndarray, optional
            A array of lambda to search; otherwise, it will be 10 ** np.linspace(-1, -7, 7).

        Returns:
        --------
        lambda_list: np.ndarray
            A array of lambda searched.

        P_list: np.ndarray
            The value defined by [Bulat, Heaton-Burgess, Cohen, Yang 2007] eqn (21).
            Corresponding to lambda in lambda_list.

        Ts_list: np.ndarray
            The Ts value for each lambda.


        """

        Ts_list = []
        L_list = []
        v_norm_list = []


        if lambda_list is None:
            lambda_list = 10 ** np.linspace(-3, -9, 7)

        if opt is None:
            opt = {
                "maxiter": 100,
                "disp": False
            }

        self.generate_components(guide_potential_components)
        self.lambda_reg = None
        # Initial calculation with no regularization
        if opt_method.lower() == 'bfgs' or opt_method.lower() == 'l-bfgs-b':
            initial_result = minimize(fun=self.lagrangian,
                                   x0=self.v0,
                                   jac=self.gradient,
                                   method=opt_method,
                                   tol=tol,
                                   options=opt
                                   )
        else:
            initial_result = minimize(fun=self.lagrangian,
                                   x0=self.v0,
                                   jac=self.gradient,
                                   hess=self.hessian,
                                   method=opt_method,
                                   tol=tol,
                                   options=opt
                                   )
        if initial_result.success == False:
            raise ValueError("Optimization was unsucessful (|grad|=%.2e) within %i iterations, "
                             "try a different intitial guess"% (np.linalg.norm(initial_result.jac), initial_result.nit)
                             + initial_result.message)
        else:
            L0 = -initial_result.fun
            initial_v0 = initial_result.x  # This is used as the initial guess for with regularization calculation.

        for reg in lambda_list:
            self.lambda_reg = reg

            if opt_method.lower() == 'bfgs' or opt_method.lower() == 'l-bfgs-b':
                opt_results = minimize(fun=self.lagrangian,
                                       x0=initial_v0,
                                       jac=self.gradient,
                                       method=opt_method,
                                       tol=tol,
                                       options=opt
                                       )

            else:
                opt_results = minimize(fun=self.lagrangian,
                                       x0=initial_v0,
                                       jac=self.gradient,
                                       hess=self.hessian,
                                       method=opt_method,
                                       tol=tol,
                                       options=opt
                                       )


            Ts_list.append(np.sum(self.T * (self.Da + self.Db)))
            v_norm_list.append(self.regul_norm)
            L_list.append(-opt_results.fun + self.lambda_reg * self.regul_norm)

        P_list = lambda_list * np.array(v_norm_list) / (L0 - np.array(L_list))

        return lambda_list, P_list, np.array(Ts_list)

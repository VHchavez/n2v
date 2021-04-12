#%%
import sys 
import time

try:    # Check for local build
    import MDI_Library as mdi
except: # Check for installed package
    import mdi


try:
    import numpy as np
    use_numpy = True
except:
    use_numpy = False

# Check for a -nompi argument
# This argument prevents the code form importing MPI
nompi_flag = False
for arg in sys.argv:
    if arg == "-nompi":
        nompi_flag = True

use_mpi4py = False
if not nompi_flag:
    try:
        from pi4py import MPI
        use_mpi4py = True
    except ImportError:
        pass

#%%

def execute_command(command, comm, self):

    if command == "EXIT":
        self.exit_flag = True

    #Bueno supongo que es aqui donde uno tiene que comunicar con el Driver
    #Entonces todo esto tiene que estar dentro del driver (????)

    # Engine <-> Driver
    # Send Commands
    elif command == "<NATOMS":
        mdi.MDI_Send(self.natoms, 1, mdi.MDI_INIT, comm)
    elif command == "<COORDS":
        mdi.MDI_Send(self.coords, 3 * self.natoms, mdi.MDI_DOUBLE, comm)

    # Receive Commands
        elif command == ">NDENSITY":
            pass
        elif command == ">CDENSITY":
            pass
        elif command == ">DENSITY":
            pass
        elif command == ">NPOTENTIAL":
            pass
        elif command == ">CPOTENTIAL":
            pass
        elif command == ">POTENTIAL":
            pass

    # Create numpy byte array?
    elif command == "<FORCES_B":
        double_size = np.dtype(np.float64).itemsize
        force_bytes = self.forces.tobytes()

        mdi.MDI_SEND(forces_bytes, 3 * self.natoms * double_size, mdi.MDI_BYTE, comm)

    else:
        raise Exception("Error in engine_py.py: MDI command not recognized")

    return 0

#%%

class MDIEngine:

    def __init__(self):
        self.exit_flag = False

        #MPI variables
        self.mpi_world = None
        self.world_rank = 0

        # Set dummy molecular information
        self.natoms = 10
        self.coords = [ 0.1 * i for i in range( 3 * self.natoms ) ]
        forces = [ 0.01 * i for i in range( 3 * self.natoms ) ]
        if use_numpy:
            self.forces = np.array(forces)
        else:
            self.forces = forces


    def run(self, mdi_options):
        # Get the MPI communicator
        if use_mpi4py:
            self.mpi_world = MPI.COMM_WORLD

        # Initialize the MDI Library
        mdi.MDI_Init(mdi_options, self.mpi_world)
        if use_mpi4py:
            self.mpi_world = mdi.MDI_MPI_get_world_comm()
            self.world_rank = self.mpi_world.Get_rank()

        # Confirm that this code is being used as an engine
        role = mdi.MDI_Get_Role()
        if not role == mdi.MDI_Engine:
            raise Exception("Must run engine.py as ENGINE")

        # Register the supported commands
        mdi.MDI_Register_None("@DEFAULT")

        # Engine '<'/'>' Driver

        # I as an engine send to Driver
        mdi.MDI_Register_Command("@DEFAUTL", "<NATOMS") 
        mdi.MDI_Register_Command("@DEFAUTL", "<COORDS") 
        mdi.MDI_Register_Command("@DEFAUTL", "<FORCES") 

        # Number of points used to represent grid. MDI_INT/1
        mdi.MDI_Register_Command("@DEFAUTL", "<NPOTENTIAL") 
        # Cartesian coordinates of a set of grid points. MDI_DOUBLE/3*NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", "<CPOTENTIAL") 
        # Set of Values of Potential on the grid. MDI_DOUBLE/NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", "<POTENTIAL")
        # Number of points used to represent grid. MDI_INT/1
        mdi.MDI_Register_Command("@DEFAUTL", "<NDENSITY") 
        # Cartesian coordinates of a set of grid points. MDI_DOUBLE/3*NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", "<CDENSITY") 
        # Set of Values of Potential on the grid. MDI_DOUBLE/NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", "<DENSITY")

        #I as an engine receive from driver
        # Number of points used to represent grid. MDI_INT/1
        mdi.MDI_Register_Command("@DEFAUTL", ">NPOTENTIAL") 
        # Cartesian coordinates of a set of grid points. MDI_DOUBLE/3*NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", ">CPOTENTIAL") 
        # Set of Values of Potential on the grid. MDI_DOUBLE/NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", ">POTENTIAL")
        # Number of points used to represent grid. MDI_INT/1
        mdi.MDI_Register_Command("@DEFAUTL", ">NDENSITY") 
        # Cartesian coordinates of a set of grid points. MDI_DOUBLE/3*NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", ">CDENSITY") 
        # Set of Values of Potential on the grid. MDI_DOUBLE/NPOTENTIAL
        mdi.MDI_Register_Command("@DEFAUTL", ">DENSITY")

        # Set the generic excecute_command function
        mdi.MDI_Set_Execute_Command_Func(execute_command, self)

        # Connect to the driver
        comm = mdi.MDI_Accept_Communicator()

        while not self.exit_flag:
            command = mdi.MDI_Recv_Command(comm)
            if use_mpi4py:
                command = self.mpi_world.bcast(command, root=0)

            execute_command( command, comm, self )

# %%

def MDI_Plugin_init_engine_py():
    engine = MDIEngine()
    engine.run("-role ENGINE -method LINK -name QM -driver_name driver")

if __name__ == "__main__":
    engine = MDIEngine()
    engine.run(sys.argv[2])
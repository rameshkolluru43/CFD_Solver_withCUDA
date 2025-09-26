import time
import sys
import CFD_Solver as CFDS
import scipy.linalg as la
import numpy as np
from scipy.sparse.linalg import spsolve,lsqr
from scipy.sparse import lil_matrix, coo_matrix
#import four_unitary_cvqls as qiskit_cvqls
# import ud4_mpi_qulacs_opp2 as qulacs_cvqls
#import ud4_qulacs_oop as qulacs_cvqls
import copy

# Example usage of compute_flux_jacobian

CFDS.readJSON("../Euler_Solver/Solver_Config.json")

# Debug Print Before Calling
print("Type of CFDS.Test_Case:", type(CFDS.Test_Case))
print("Test case read from JSON:", CFDS.Test_Case)  

# Ensure it's an integer
if isinstance(CFDS.Test_Case, int):
    testcase = CFDS.Test_Case
else:
    print("Error: CFDS.Test_Case is not an integer!")
    testcase = int(CFDS.Test_Case)  # Try to convert if needed

# Debug Print
print("Passing test case:", testcase)

# Call the C++ function
CFDS.testCase(testcase)

iterations = 0
Total_Iterations = 1000

interval = Total_Iterations       # Frequency at which CVQLS Solver is used. If it is set to Total_Iterations then, last iteration is run on CVQLS
CVQLS_Implementation = "Qlacs"
Fine_Tune = False
optimizer = "COBYLA"

# Main loop
nCells = CFDS.get_NoPhysical_Cells()
    
print("ncells =",nCells)


Grid_Vtk_File = CFDS.Get_Grid_VtkFile()
Error_File = CFDS.Get_ErrorFileName()
Initial_Solution_File =CFDS.Get_Initial_Solution_FileName()
Solution_File = CFDS.Get_SolutionFile()
Final_Solution_File =CFDS.Get_Final_Solution_FileName()

print(Grid_Vtk_File)
print(Solution_File)
print(Final_Solution_File)

CFDS.A= [[0.0 for _ in range(4*nCells)] for _ in range(4*nCells)]  # Initialize a 4x4 matrix
CFDS.b= [0.0 for _ in range(4*nCells)]  # Initialize a 4x4 matrix
A = lil_matrix((4 * nCells, 4 * nCells))

# Assuming you have the following three lists returned from the C++ side
total_size = 4*nCells     # The size of the matrix (4 * Total_No_Cells)

Solution_Data_Type = 1
for iterations in range(1, Total_Iterations + 1):
    # timer equivalent of clock in Python
#    timer_start = time.perf_counter()

    # Apply boundary conditions
    CFDS.Apply_Boundary_Conditions()
    
    b = CFDS.Assemble_b(CFDS.b)
    #print(len(b))

    #print("Assembling b done")
    
    min_dt = CFDS.get_Min_dt()
    #print(min_dt)
    
    if iterations%interval == 0:
        print("Entered into CVQLS Solver")
        A_full =  CFDS.Assemble_A(CFDS.A, min_dt);
        if CVQLS_Implementation == "Qiskit":
            SolObj = qiskit_cvqls.Qsolver(np.array(A_full),np.array(b))
            x = SolObj.solve(optimizer,Fine_Tune)
        elif CVQLS_Implementation == "Qulacs":
            SolObj = qulacs_cvqls.Qsolver(np.array(A_full),np.array(b))
            x = SolObj.solve(optimizer,Fine_Tune)
    else:
        CFDS.Assemble_A1(min_dt)
        row_indices = CFDS.get_row_indices()  # List of row indices
        
        col_indices = CFDS.get_col_indices()  # List of column indices
        
        values = CFDS.get_Values()       # List of non-zero values
        A_coo = coo_matrix((values, (row_indices, col_indices)), shape=(4*nCells, 4*nCells))
        #  Convert to CSR matrix
        A_csr = A_coo.tocsr()
        x = lsqr(A_csr, b)[0]
        x_prev = copy.copy(x)
        
    #print("Solving Ax = b done")
    CFDS.Set_DelU(x)
    #print("setting delU to c++ code")

    # Estimate error and update solution
    CFDS.Estimate_Error()
    #print("estimating the error")
    CFDS.Update()
    #print("updating the solution vector")
    # Write solution every 100 iterations
    CFDS.iterations=iterations
    #timer_start = time.perf_counter()
    if iterations % 100 == 0:
        CFDS.WriteErrorFile(Error_File)
        CFDS.Write_Solution_1(Solution_File, Solution_Data_Type)
        CFDS.Read_Write_Grid(Grid_Vtk_File, Final_Solution_File)
        CFDS.Append_Solution(Solution_File, Final_Solution_File)
        print(iterations)

    #end_time = time.time()

quantum_error = x_prev - x
print("Difference between the Quantum and Classical results = ", max(abs(quantum_error)))

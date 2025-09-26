import CFD_Solver as CFDS
import scipy.linalg as la

# Example usage of compute_flux_jacobian

CFDS.SetGlobalVariables(
        grid_type=0,
        init_type=0,
        limiter_case=0,
        area_weighted_avg=1,
        flux_type=1,
        is_second_order=True,
        time_accurate=False,
        local_time_stepping=False,
        non_dimensional_form=False,
        is_weno=False,
        is_char=False,
        dissipation_type=1,
        is_movers_1=False,
        enable_entropy_fix=False,
        cfl=0.1,
        test_case=3,  # Ramp 15 Degree case as an example
        grid_size=2,
        iterations=100000,
        Implicit_Method=True
    )
CFDS.Ramp_15_Degree()  # Call the function

CFDS.Cell_No = 1
CFDS.Face_No = 2
CFDS.Ac = [[0.0 for _ in range(4)] for _ in range(4)]  # Initialize a 4x4 matrix
print("Computed Flux Jacobian:", CFDS.Ac)

updated_Ac = CFDS.Compute_Flux_Jacobian(CFDS.Cell_No, CFDS.Ac, CFDS.Face_No)

print("Computed Flux Jacobian:", updated_Ac)
# Print the matrix to check if it's updated
print("Updated Ac matrix:")
for row in updated_Ac:
    print(row)
    
nCells = CFDS.get_NoPhysical_Cells()
    
print(nCells)
    
CFDS.A= [[0.0 for _ in range(4*nCells)] for _ in range(4*nCells)]  # Initialize a 4x4 matrix
A = CFDS.Assemble_A(CFDS.A,0.01)

#for row in A:
#    print(row)
CFDS.b= [0.0 for _ in range(4*nCells)]  # Initialize a 4x4 matrix
b = CFDS.Assemble_b(CFDS.b)

#print(b)

x = la.solve(A,b)

print("Solution Vector")
#print(x)

CFDS.Set_DelU(x)

CFDS.Estimate_Error()
CFDS.Update()
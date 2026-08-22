import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from scipy.sparse import isspmatrix

def print_matrix_pattern(matrix, filled_symbol="■", zero_symbol="."):
    """
    Prints a pattern representation of a matrix:
    - filled_symbol for non-zero entries
    - zero_symbol for zero entries

    Parameters:
        matrix (np.ndarray or scipy.sparse matrix): The matrix to visualize.
        filled_symbol (str): The symbol for non-zero entries (default: "■").
        zero_symbol (str): The symbol for zero entries (default: ".").
    """
    # Convert sparse matrix to dense if needed
    if isspmatrix(matrix):
        matrix = matrix.toarray()  # Convert to NumPy array

    if not isinstance(matrix, np.ndarray):
        raise ValueError("The input matrix must be a NumPy array or a sparse matrix.")

    rows, cols = matrix.shape
    print(f"Matrix shape: {rows}x{cols}")
    
    for row in matrix:
        line = "".join(filled_symbol if val != 0 else zero_symbol for val in row)
        print(line)

def plot_matrix_pattern(matrix, title="Matrix Sparsity Pattern"):
    """
    Plots the sparsity pattern of a given matrix.

    Parameters:
        matrix (scipy.sparse or numpy.ndarray): The matrix to plot.
        title (str): Title of the plot.
    """
    # Convert to dense if the input is a sparse matrix
    if not isinstance(matrix, np.ndarray):
        matrix = matrix.toarray()

    plt.figure(figsize=(8, 8))
    plt.spy(matrix, markersize=1)
    plt.title(title, fontsize=16)
    plt.xlabel("Column Index", fontsize=12)
    plt.ylabel("Row Index", fontsize=12)
    plt.grid(False)
    plt.show()


def f(x, y, is_poisson=True):
    """Source term for Poisson equation. Returns zero for Laplace."""
    if is_poisson:
        return np.sin(np.pi * x) * np.sin(np.pi * y)
    else:
        return np.zeros_like(x)

def apply_dirichlet_boundary(A, b, nx, ny, dirichlet_values):
    """Apply Dirichlet boundary conditions."""
    for i in range(nx):
        # Bottom boundary
        A[i, :] = 0
        A[i, i] = 1
        b[i] = dirichlet_values['bottom']
        
        # Top boundary
        top_index = (ny - 1) * nx + i
        A[top_index, :] = 0
        A[top_index, top_index] = 1
        b[top_index] = dirichlet_values['top']
    
    for j in range(ny):
        # Left boundary
        left_index = j * nx
        A[left_index, :] = 0
        A[left_index, left_index] = 1
        b[left_index] = dirichlet_values['left']
        
        # Right boundary
        right_index = j * nx + (nx - 1)
        A[right_index, :] = 0
        A[right_index, right_index] = 1
        b[right_index] = dirichlet_values['right']

def apply_neumann_boundary(b, nx, ny, dx, dy, neumann_gradients):
    """Apply Neumann boundary conditions."""
    for i in range(nx):
        # Top boundary
        top_index = (ny - 1) * nx + i
        b[top_index] -= neumann_gradients['top'] * dy
        
        # Bottom boundary
        b[i] -= neumann_gradients['bottom'] * dy
    
    for j in range(ny):
        # Left boundary
        left_index = j * nx
        b[left_index] -= neumann_gradients['left'] * dx
        
        # Right boundary
        right_index = j * nx + (nx - 1)
        b[right_index] -= neumann_gradients['right'] * dx

def solve_fdm(nx, ny, BC_Method, is_poisson=True):
    """Solve the Poisson or Laplace equation using FDM."""
    lx, ly = 1.0, 1.0
    dx, dy = lx / (nx - 1), ly / (ny - 1)
    x = np.linspace(0, lx, nx)
    y = np.linspace(0, ly, ny)
    xx, yy = np.meshgrid(x, y)
    f_values = f(xx, yy, is_poisson).flatten()

    n = nx * ny

    # Sparse matrix construction
    main_diag = -4 * np.ones(n)
    side_diag = np.ones(n - 1)
    side_diag[np.arange(1, n) % nx == 0] = 0  # Boundary adjustment
    up_down_diag = np.ones(n - nx)
    diagonals = [main_diag, side_diag, side_diag, up_down_diag, up_down_diag]
    offsets = [0, -1, 1, -nx, nx]
    A = diags(diagonals, offsets, shape=(n, n)).tocsc()
    
    # RHS vector
    b = -f_values

    # Apply Boundary Conditions
    if BC_Method == 1:
        # Dirichlet Boundary Conditions
        dirichlet_values = {'top': 1, 'bottom': 1, 'left': 1, 'right': 1}
        apply_dirichlet_boundary(A, b, nx, ny, dirichlet_values)
    elif BC_Method == 2:
        # Neumann Boundary Conditions
        neumann_gradients = {'top': 1, 'bottom': -1, 'left': 0, 'right': 0}
        apply_neumann_boundary(b, nx, ny, dx, dy, neumann_gradients)
    elif BC_Method == 3:
        # Mixed Boundary Conditions
        dirichlet_values = {'top': 1, 'bottom': 0, 'left': 0, 'right': 1}
        neumann_gradients = {'top': 0, 'bottom': 0, 'left': 0, 'right': 0}
        apply_dirichlet_bondary(A, b, nx, ny, dirichlet_values)
        apply_neumann_boundary(b, nx, ny, dx, dy, neumann_gradients)
        
        
    #print_matrix_pattern(A, filled_symbol="■", zero_symbol=".")
    #plot_matrix_pattern(A, title="Matrix Sparsity Pattern")
   
    
    # Solve system
    u = spsolve(A, b).reshape((ny, nx))
    equation_type = "Poisson" if is_poisson else "Laplace"
    plot_solution_xy(xx, yy, u, f"{equation_type} Equation Solution (FDM)", is_triangulation=False)

def plot_solution_xy(x, y, u, title, is_triangulation=False):
    """Plot the solution of the equation."""
    plt.figure(figsize=(8, 6))
    plt.contourf(x, y, u, levels=50, cmap='viridis')
    plt.colorbar(label='Solution u(x, y)')
    plt.title(title)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()

# Example Usage
if __name__ == "__main__":
    nx, ny = 5, 5
    BC_Method = 1  # 1: Dirichlet, 2: Neumann, 3: Mixed
    is_poisson = True  # Change to True for Poisson equation
    solve_fdm(nx, ny, BC_Method, is_poisson)
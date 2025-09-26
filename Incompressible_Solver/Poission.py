import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.spatial import Delaunay

# Source term for Poisson's equation
f = lambda x, y: np.sin(np.pi * x) * np.sin(np.pi * y)

def plot_solution_xy(x, y, u, method_name, is_triangulation=False, triangles=None):
    """
    Plot the solution of the Poisson equation.

    Args:
        x (ndarray): 2D array of x-coordinates for structured grids or 1D for unstructured.
        y (ndarray): 2D array of y-coordinates for structured grids or 1D for unstructured.
        u (ndarray): Solution values (2D for structured, 1D for unstructured).
        method_name (str): Name of the method for the title.
        is_triangulation (bool): True if the solution uses a triangulated grid.
        triangles (ndarray): Connectivity for triangulated grids.
    """
    plt.figure(figsize=(6, 5))
    if is_triangulation:
        plt.tricontourf(x, y, triangles, u, levels=50, cmap='viridis')
    else:
        plt.contourf(x, y, u, levels=50, cmap='viridis')  # FDM/structured grids
    plt.colorbar(label='u')
    plt.title(f'{method_name} Solution to Poisson Equation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    
def plot_solution(points, grid_shape, u, method_name, is_triangulation=False, triangles=None):
    plt.figure(figsize=(6, 5))
    if is_triangulation:
        plt.tricontourf(points[:, 0], points[:, 1], triangles, u, levels=50, cmap='viridis')
    else:
        # Reshape the solution for quadrilateral grid visualization
        x = points[:, 0].reshape(grid_shape)
        y = points[:, 1].reshape(grid_shape)
        u = u.reshape(grid_shape)
        plt.contourf(x, y, u, levels=50, cmap='viridis')
    plt.colorbar(label='u')
    plt.title(f'{method_name} Solution to Poisson Equation')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()    
    
def apply_dirichlet_boundary(grid, boundary_values):
    """
    Applies Dirichlet boundary conditions to the grid.
    
    Parameters:
    grid : 2D numpy array
        The grid representing the solution.
    boundary_values : dict
        Dictionary with keys ['top', 'bottom', 'left', 'right'] specifying
        the Dirichlet values for the respective boundaries.
    """
    nx, ny = grid.shape
    
    # Top boundary (y = ny-1)
    grid[0, :] = boundary_values['top']
    
    # Bottom boundary (y = 0)
    grid[-1, :] = boundary_values['bottom']
    
    # Left boundary (x = 0)
    grid[:, 0] = boundary_values['left']
    
    # Right boundary (x = nx-1)
    grid[:, -1] = boundary_values['right']

def apply_neumann_boundary(grid, dx, dy, boundary_gradients):
    """
    Applies Neumann boundary conditions to the grid.
    
    Parameters:
    grid : 2D numpy array
        The grid representing the solution.
    dx, dy : float
        The grid spacing in the x and y directions.
    boundary_gradients : dict
        Dictionary with keys ['top', 'bottom', 'left', 'right'] specifying
        the Neumann gradient values (du/dn) for the respective boundaries.
    """
    nx, ny = grid.shape
    
    # Top boundary (y = ny-1)
    grid[0, 1:-1] = grid[1, 1:-1] + boundary_gradients['top'] * dy
    
    # Bottom boundary (y = 0)
    grid[-1, 1:-1] = grid[-2, 1:-1] - boundary_gradients['bottom'] * dy
    
    # Left boundary (x = 0)
    grid[1:-1, 0] = grid[1:-1, 1] - boundary_gradients['left'] * dx
    
    # Right boundary (x = nx-1)
    grid[1:-1, -1] = grid[1:-1, -2] + boundary_gradients['right'] * dx

def apply_mixed_boundary(grid, dx, dy, dirichlet_values, neumann_gradients):
    """
    Applies mixed boundary conditions to the grid.
    
    Parameters:
    grid : 2D numpy array
        The grid representing the solution.
    dx, dy : float
        The grid spacing in the x and y directions.
    dirichlet_values : dict
        Dictionary with keys ['top', 'bottom', 'left', 'right'] specifying
        the Dirichlet values for respective boundaries.
    neumann_gradients : dict
        Dictionary with keys ['top', 'bottom', 'left', 'right'] specifying
        the Neumann gradient values for respective boundaries.
    """
    nx, ny = grid.shape
    
    # Top boundary (mixed: Dirichlet on edges, Neumann otherwise)
    grid[0, 1:-1] = grid[1, 1:-1] + neumann_gradients['top'] * dy
    grid[0, 0] = dirichlet_values['top']  # Top-left corner
    grid[0, -1] = dirichlet_values['top']  # Top-right corner
    
    # Bottom boundary (mixed)
    grid[-1, 1:-1] = grid[-2, 1:-1] - neumann_gradients['bottom'] * dy
    grid[-1, 0] = dirichlet_values['bottom']  # Bottom-left corner
    grid[-1, -1] = dirichlet_values['bottom']  # Bottom-right corner
    
    # Left boundary (mixed)
    grid[1:-1, 0] = grid[1:-1, 1] - neumann_gradients['left'] * dx
    grid[0, 0] = dirichlet_values['left']  # Top-left corner
    grid[-1, 0] = dirichlet_values['left']  # Bottom-left corner
    
    # Right boundary (mixed)
    grid[1:-1, -1] = grid[1:-1, -2] + neumann_gradients['right'] * dx
    grid[0, -1] = dirichlet_values['right']  # Top-right corner
    grid[-1, -1] = dirichlet_values['right']  # Bottom-right corner

def solve_poisson_fdm(nx, ny, BC_Method):
    lx, ly = 1.0, 1.0
    dx, dy = lx / (nx - 1), ly / (ny - 1)
    x = np.linspace(0, lx, nx)
    y = np.linspace(0, ly, ny)
    xx, yy = np.meshgrid(x, y)
    f_values = f(xx, yy).flatten()

    n = nx * ny

    # Sparse matrix construction
    main_diag = -4 * np.ones(n)
    side_diag = np.ones(n - 1)
    side_diag[np.arange(1, n) % nx == 0] = 0  # Boundary adjustment
    up_down_diag = np.ones(n - nx)
    diagonals = [main_diag, side_diag, side_diag, up_down_diag, up_down_diag]
    offsets = [0, -1, 1, -nx, nx]
    A = diags(diagonals, offsets, shape=(n, n)).tocsc()

    # RHS and boundary conditions
    b = -f_values
    
    # Apply Boundary conditions 
    if BC_Method == 1:
        # Dirichlet boundary values
        dirichlet_values = {'top': 1, 'bottom': 0, 'left': 0, 'right': 1}
        apply_dirichlet_boundary(A, dirichlet_values)
        print("Grid with Dirichlet Boundary Conditions:")
        print(grid)
        
    elif BC_Method ==2:
        # Reset grid and apply Neumann conditions
            grid = np.zeros((nx, ny))
            neumann_gradients = {'top': 1, 'bottom': -1, 'left': 0, 'right': 0}
            apply_neumann_boundary(grid, dx=0.1, dy=0.1, boundary_gradients=neumann_gradients)
            print("\nGrid with Neumann Boundary Conditions:")
            print(grid)
    elif BC_Method == 3 :
        # Mixed boundary conditions
            grid = np.zeros((nx, ny))
            apply_mixed_boundary(grid, dx=0.1, dy=0.1, 
            dirichlet_values={'top': 1, 'bottom': 0, 'left': 0, 'right': 1},
            neumann_gradients={'top': 1, 'bottom': -1, 'left': 0, 'right': 0})
            print("\nGrid with Mixed Boundary Conditions:")
            print(grid)
    
    for i in range(nx):
        b[i] = 0  # Bottom boundary
        b[(ny - 1) * nx + i] = 0  # Top boundary
    for j in range(ny):
        b[j * nx] = 0  # Left boundary
        b[j * nx + nx - 1] = 0  # Right boundary

    # Solve system
    u = spsolve(A, b).reshape((ny, nx))
    plot_solution_xy(xx, yy, u, "FDM", is_triangulation=False)

# Finite Element Method (FEM)
def solve_poisson_fem(points, triangles):
    n_nodes = len(points)
    A = np.zeros((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    for tri in triangles:
        p1, p2, p3 = points[tri]
        area = 0.5 * abs(np.linalg.det([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]]))
        grads = np.linalg.inv([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]])[1:]

        K_local = area * np.dot(grads.T, grads)

        for i in range(3):
            for j in range(3):
                A[tri[i], tri[j]] += K_local[i, j]

        centroid = np.mean([p1, p2, p3], axis=0)
        b_local = area / 3 * f(centroid[0], centroid[1])
        for i in range(3):
            b[tri[i]] += b_local

    boundary_nodes = [i for i, (x, y) in enumerate(points) if x == 0 or x == 1 or y == 0 or y == 1]
    A[boundary_nodes, :] = 0
    A[boundary_nodes, boundary_nodes] = 1
    b[boundary_nodes] = 0

    u = np.linalg.solve(A, b)
    plot_solution(points, None, u, "FEM", is_triangulation=True, triangles=triangles)

def solve_poisson_fem_quads(points,quads):
    """
    Solve the 2D Poisson equation using FEM on a quadrilateral mesh.

    Args:
        points (ndarray): Array of node coordinates (N_nodes, 2).
        quads (ndarray): Array of quadrilateral elements (N_elements, 4) as indices into the points array.
    """
    n_nodes = len(points)
    A = np.zeros((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    # Define Gaussian quadrature points and weights for a quadrilateral
    gauss_points = [
        (-np.sqrt(1/3), -np.sqrt(1/3)),
        ( np.sqrt(1/3), -np.sqrt(1/3)),
        ( np.sqrt(1/3),  np.sqrt(1/3)),
        (-np.sqrt(1/3),  np.sqrt(1/3)),
    ]
    weights = [1, 1, 1, 1]

    # Loop over each quadrilateral
    for quad in quads:
        p1, p2, p3, p4 = points[quad]
        vertices = np.array([p1, p2, p3, p4])

        # Jacobian matrix and determinant for mapping to physical space
        def jacobian(xi, eta):
            dN_dxi = np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)]]) / 4
            dN_deta = np.array([[-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]]) / 4
            J = np.dot(np.vstack((dN_dxi, dN_deta)), vertices)
            return J, np.linalg.det(J)

        # Shape functions for the bilinear quadrilateral
        def shape_functions(xi, eta):
            return np.array([
                (1 - xi) * (1 - eta) / 4,
                (1 + xi) * (1 - eta) / 4,
                (1 + xi) * (1 + eta) / 4,
                (1 - xi) * (1 + eta) / 4,
            ])

        # Loop over Gauss points
        K_local = np.zeros((4, 4))
        b_local = np.zeros(4)
        for gp, w in zip(gauss_points, weights):
            xi, eta = gp
            J, detJ = jacobian(xi, eta)
            invJ = np.linalg.inv(J)

            # Compute gradients of shape functions in physical space
            dN_dxi = np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)]]) / 4
            dN_deta = np.array([[-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]]) / 4
            dN_dx_dy = np.dot(invJ, np.vstack((dN_dxi, dN_deta)))

            # Local stiffness matrix
            K_local += w * detJ * np.dot(dN_dx_dy.T, dN_dx_dy)

            # Local load vector
            N = shape_functions(xi, eta)
            x, y = np.dot(N, vertices)  # Map Gauss point to physical space
            b_local += w * detJ * N * f(x, y)

        # Assemble local contributions into the global matrix and vector
        for i in range(4):
            for j in range(4):
                A[quad[i], quad[j]] += K_local[i, j]
            b[quad[i]] += b_local[i]

    # Apply Dirichlet boundary conditions (assuming u = 0 on the boundary)
    boundary_nodes = [i for i, (x, y) in enumerate(points) if x == 0 or x == 1 or y == 0 or y == 1]
    A[boundary_nodes, :] = 0
    A[boundary_nodes, boundary_nodes] = 1
    b[boundary_nodes] = 0

    # Solve the linear system
    u = np.linalg.solve(A, b)

    # Plot the solution
    ny = int(np.sqrt(len(points)))  # Number of rows in the grid
    nx = len(points) // ny          # Number of columns in the grid
    plot_solution(points, (ny, nx), u, "FEM (Quadrilateral Mesh)", is_triangulation=False)
    
# Finite Volume Method (FVM)
def solve_poisson_fvm(nx, ny):
    lx, ly = 1.0, 1.0
    dx, dy = lx / (nx - 1), ly / (ny - 1)
    x = np.linspace(0, lx, nx)
    y = np.linspace(0, ly, ny)
    xx, yy = np.meshgrid(x, y)
    f_values = f(xx, yy).flatten()

    n = nx * ny

    main_diag = -4 * np.ones(n)
    side_diag = np.ones(n - 1)
    side_diag[np.arange(1, n) % nx == 0] = 0
    up_down_diag = np.ones(n - nx)
    diagonals = [main_diag, side_diag, side_diag, up_down_diag, up_down_diag]
    offsets = [0, -1, 1, -nx, nx]
    A = diags(diagonals, offsets, shape=(n, n)).tocsc()

    b = -f_values
    for i in range(nx):
        b[i] = 0
        b[(ny - 1) * nx + i] = 0
    for j in range(ny):
        b[j * nx] = 0
        b[j * nx + nx - 1] = 0

    u = spsolve(A, b).reshape((ny, nx))
    plot_solution(xx, yy, u, "FVM")

def solve_poisson_fvm_unstructured(points, triangles):
    """
    Solve the 2D Poisson equation using FVM on arbitrary-oriented control volumes.

    Args:
        points (ndarray): Array of node coordinates (N_nodes, 2).
        triangles (ndarray): Array of triangular elements (N_elements, 3) as indices into the points array.
    """
    n_nodes = len(points)
    n_elements = len(triangles)

    # Initialize the stiffness matrix and load vector
    A = np.zeros((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    # Loop over each triangle to assemble contributions
    for tri in triangles:
        # Get vertex coordinates
        p1, p2, p3 = points[tri]
        vertices = np.array([p1, p2, p3])

        # Compute the area of the triangle
        area = 0.5 * abs(np.linalg.det([[1, p1[0], p1[1]], [1, p2[0], p2[1]], [1, p3[0], p3[1]]]))

        # Compute edge vectors and outward normals
        edges = np.array([p2 - p1, p3 - p2, p1 - p3])
        normals = np.array([[-edge[1], edge[0]] for edge in edges])  # Rotate edges by 90 degrees
        normals /= np.linalg.norm(normals, axis=1)[:, np.newaxis]  # Normalize

        # Local stiffness contributions (fluxes)
        for i, edge_normal in enumerate(normals):
            v_start, v_end = vertices[i], vertices[(i + 1) % 3]
            edge_length = np.linalg.norm(v_end - v_start)

            # Flux contribution for edge
            for a, b in zip([i, (i + 1) % 3], [(i + 1) % 3, i]):
                contribution = np.dot(edge_normal, (v_end + v_start) / 2 - points[tri[a]]) * edge_length
                A[tri[a], tri[b]] += contribution / 2
                A[tri[b], tri[a]] += contribution / 2

        # Source term contribution
        centroid = np.mean(vertices, axis=0)
        f_centroid = f(centroid[0], centroid[1])  # Evaluate source term at centroid
        for i in range(3):
            b[tri[i]] += f_centroid * area / 3

    # Apply Dirichlet boundary conditions (assume u = 0 on the boundary)
    boundary_nodes = [i for i, (x, y) in enumerate(points) if x == 0 or x == 1 or y == 0 or y == 1]
    A[boundary_nodes, :] = 0
    A[boundary_nodes, boundary_nodes] = 1
    b[boundary_nodes] = 0

    # Solve the linear system
    u = np.linalg.solve(A, b)

    # Plot the solution
    plot_solution(points, None, u, "FVM (Unstructured)", is_triangulation=True, triangles=triangles)
    
def solve_poisson_fvm_quad_mesh(nx, ny):
    """
    Generate a quadrilateral mesh and solve the 2D Poisson equation using FVM.

    Args:
        nx (int): Number of points in the x-direction.
        ny (int): Number of points in the y-direction.
    """
    # Step 1: Generate the mesh points
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    grid_x, grid_y = np.meshgrid(x, y)
    points = np.c_[grid_x.ravel(), grid_y.ravel()]

    # Step 2: Create quadrilateral connectivity
    quads = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            p1 = j * nx + i
            p2 = p1 + 1
            p3 = p2 + nx
            p4 = p1 + nx
            quads.append([p1, p2, p3, p4])
    quads = np.array(quads)

    # Step 3: Solve the Poisson equation using FVM
    solve_poisson_fvm_quads(points, quads)
        
def solve_poisson_fvm_tri_mesh(nx, ny):
    """
    Generate a triangular mesh and solve the 2D Poisson equation using FVM.

    Args:
        nx (int): Number of points in the x-direction.
        ny (int): Number of points in the y-direction.
    """
    # Step 1: Generate the mesh points
    x = np.linspace(0, 1, nx)
    y = np.linspace(0, 1, ny)
    grid_x, grid_y = np.meshgrid(x, y)
    points = np.c_[grid_x.ravel(), grid_y.ravel()]

    # Step 2: Create a triangular mesh using Delaunay triangulation
    tri = Delaunay(points)

    # Step 3: Solve the Poisson equation using FVM
    solve_poisson_fvm_unstructured(points, tri.simplices)

def solve_poisson_fvm_quads(points, quads):
    """
    Solve the 2D Poisson equation using FVM on a quadrilateral mesh.
    
    Args:
        points (ndarray): Array of node coordinates (N_nodes, 2).
        quads (ndarray): Array of quadrilateral elements (N_elements, 4) as indices into the points array.
    """
    # Similar to solve_poisson_fvm_unstructured but adapted for quadrilaterals
    n_nodes = len(points)
    A = np.zeros((n_nodes, n_nodes))
    b = np.zeros(n_nodes)

    for quad in quads:
        p1, p2, p3, p4 = points[quad]
        vertices = np.array([p1, p2, p3, p4])

        # Compute centroid and area (approximation for FVM)
        centroid = np.mean(vertices, axis=0)
        f_centroid = f(centroid[0], centroid[1])  # Evaluate source term at centroid
        for node in quad:
            b[node] += f_centroid / 4  # Distribute source term equally to nodes

        # Local stiffness contributions (fluxes)
        for i, node in enumerate(quad):
            next_node = quad[(i + 1) % 4]
            edge_vector = points[next_node] - points[node]
            edge_length = np.linalg.norm(edge_vector)
            normal = np.array([-edge_vector[1], edge_vector[0]]) / edge_length  # Outward normal
            A[node, node] += edge_length  # Diagonal contribution
            A[node, next_node] -= edge_length  # Off-diagonal contribution

    # Apply Dirichlet boundary conditions
    boundary_nodes = [i for i, (x, y) in enumerate(points) if x == 0 or x == 1 or y == 0 or y == 1]
    A[boundary_nodes, :] = 0
    A[boundary_nodes, boundary_nodes] = 1
    b[boundary_nodes] = 0

    # Solve the system
    u = np.linalg.solve(A, b)

    # Plot the solution
    ny = int(np.sqrt(len(points)))
    nx = len(points) // ny
    plot_solution(points, (ny, nx), u, "FVM (Quadrilateral Mesh)", is_triangulation=False)
        
            
def main():
    """
    Main function to solve the 2D Poisson equation using different methods (FDM, FEM, FVM)
    and mesh types (triangular or quadrilateral).
    """
    print("Select a method to solve the 2D Poisson Equation:")
    print("1. Finite Difference Method (FDM)")
    print("2. Finite Element Method (FEM)")
    print("3. Finite Volume Method (FVM)")
    method = input("Enter your choice (1/2/3): ")
    
    print("Select Boundary Conditions to be applied to solve the 2D Poisson Equation:")
    print("1. Dirchlet Conditions")
    print("2. Naumann Conditions")
    print("3. Mixed or Robin Conditions ")
    BC_Method = input("Enter your choice (1/2/3): ")
    
    # Prompt the user for input and process it
    Nx, Ny = map(int, input("Enter Number of Points in x and y Directions (e.g., 10 10): ").split())    
    
    if method == "1":
        # Solve using FDM
        print("\nSolving using FDM...")
        solve_poisson_fdm(Nx, Ny,BC_Method)

    elif method == "2":
        # Solve using FEM
        print("\nSelect a mesh type for FEM:")
        print("1. Triangular mesh")
        print("2. Quadrilateral mesh")
        mesh_choice = input("Enter your choice (1/2): ")

        if mesh_choice == "1":
            # Solve FEM with triangular mesh
            print("\nSolving FEM with triangular mesh...")
            nx, ny = 10, 10
            x = np.linspace(0, 1, nx)
            y = np.linspace(0, 1, ny)
            grid_x, grid_y = np.meshgrid(x, y)
            points = np.c_[grid_x.ravel(), grid_y.ravel()]
            tri = Delaunay(points)
            solve_poisson_fem(points, tri.simplices)

        elif mesh_choice == "2":
            # Solve FEM with quadrilateral mesh
            print("\nSolving FEM with quadrilateral mesh...")
            nx, ny = 10, 10
            x = np.linspace(0, 1, nx)
            y = np.linspace(0, 1, ny)
            grid_x, grid_y = np.meshgrid(x, y)
            points = np.c_[grid_x.ravel(), grid_y.ravel()]
            quads = []
            for j in range(ny - 1):
                for i in range(nx - 1):
                    p1 = j * nx + i
                    p2 = p1 + 1
                    p3 = p2 + nx
                    p4 = p1 + nx
                    quads.append([p1, p2, p3, p4])
            quads = np.array(quads)
            solve_poisson_fem_quads(points, quads)  # Use the correct function for quadrilaterals
        else:
            print("Invalid choice for FEM mesh type.")

    elif method == "3":
        # Solve using FVM
        print("\nSelect a mesh type for FVM:")
        print("1. Triangular mesh")
        print("2. Quadrilateral mesh")
        mesh_choice = input("Enter your choice (1/2): ")

        if mesh_choice == "1":
            # Solve FVM with triangular mesh
            print("\nSolving FVM with triangular mesh...")
            solve_poisson_fvm_tri_mesh(Nx, Ny)

        elif mesh_choice == "2":
            # Solve FVM with quadrilateral mesh
            print("\nSolving FVM with quadrilateral mesh...")
            solve_poisson_fvm_quad_mesh(Nx, Ny)

        else:
            print("Invalid choice for FVM mesh type.")

    else:
        print("Invalid choice for the solution method.")

# Execute main
if __name__ == "__main__":
    main()

    
    
    
    

    
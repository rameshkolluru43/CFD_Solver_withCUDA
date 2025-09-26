/**
 * @file cuda_kernels.cuh
 * @brief CUDA kernel declarations for CFD computations
 * @author Ramesh Kolluru
 * @date 2024
 * @version 0.1.0-alpha
 * 
 * This file contains the declarations for all CUDA kernels used in the
 * CFD solver. These kernels implement various numerical operations
 * optimized for GPU execution.
 */

#ifndef CUDA_KERNELS_CUH
#define CUDA_KERNELS_CUH

#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/**
 * @namespace cfd
 * @brief Main namespace for CFD solver components
 */
namespace cfd {

/**
 * @namespace kernels
 * @brief Namespace containing CUDA kernel functions
 */
namespace kernels {

/**
 * @brief CUDA kernel for computing velocity divergence
 * @param u_d Device pointer to x-velocity field
 * @param v_d Device pointer to y-velocity field
 * @param w_d Device pointer to z-velocity field (optional, can be nullptr for 2D)
 * @param div_d Device pointer to output divergence field
 * @param nx Number of grid points in x-direction
 * @param ny Number of grid points in y-direction
 * @param nz Number of grid points in z-direction
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @param dz Grid spacing in z-direction
 * 
 * This kernel computes the divergence of the velocity field using
 * finite differences. For 2D problems, set w_d to nullptr and nz to 1.
 */
__global__ void computeDivergence(const double* u_d, const double* v_d, const double* w_d,
                                 double* div_d,
                                 int nx, int ny, int nz,
                                 double dx, double dy, double dz);

/**
 * @brief CUDA kernel for computing pressure gradient
 * @param p_d Device pointer to pressure field
 * @param grad_px_d Device pointer to output x-component of pressure gradient
 * @param grad_py_d Device pointer to output y-component of pressure gradient
 * @param grad_pz_d Device pointer to output z-component of pressure gradient (optional)
 * @param nx Number of grid points in x-direction
 * @param ny Number of grid points in y-direction
 * @param nz Number of grid points in z-direction
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @param dz Grid spacing in z-direction
 * 
 * Computes the gradient of the pressure field using central differences.
 */
__global__ void computePressureGradient(const double* p_d,
                                       double* grad_px_d, double* grad_py_d, double* grad_pz_d,
                                       int nx, int ny, int nz,
                                       double dx, double dy, double dz);

/**
 * @brief CUDA kernel for computing viscous terms (Laplacian)
 * @param field_d Device pointer to input field
 * @param laplacian_d Device pointer to output Laplacian
 * @param nx Number of grid points in x-direction
 * @param ny Number of grid points in y-direction
 * @param nz Number of grid points in z-direction
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @param dz Grid spacing in z-direction
 * @param viscosity Kinematic viscosity
 * 
 * Computes the viscous terms using a second-order finite difference
 * approximation of the Laplacian operator.
 */
__global__ void computeViscousTerms(const double* field_d, double* laplacian_d,
                                   int nx, int ny, int nz,
                                   double dx, double dy, double dz,
                                   double viscosity);

/**
 * @brief CUDA kernel for convective terms computation
 * @param u_d Device pointer to x-velocity field
 * @param v_d Device pointer to y-velocity field
 * @param w_d Device pointer to z-velocity field
 * @param conv_u_d Device pointer to output convective term for u
 * @param conv_v_d Device pointer to output convective term for v
 * @param conv_w_d Device pointer to output convective term for w
 * @param nx Number of grid points in x-direction
 * @param ny Number of grid points in y-direction
 * @param nz Number of grid points in z-direction
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @param dz Grid spacing in z-direction
 * 
 * Computes the nonlinear convective terms (u∇u) using upwind differencing
 * for stability.
 */
__global__ void computeConvectiveTerms(const double* u_d, const double* v_d, const double* w_d,
                                      double* conv_u_d, double* conv_v_d, double* conv_w_d,
                                      int nx, int ny, int nz,
                                      double dx, double dy, double dz);

/**
 * @brief CUDA kernel for applying boundary conditions
 * @param field_d Device pointer to field array
 * @param boundary_type Type of boundary condition (0=Dirichlet, 1=Neumann, 2=Periodic)
 * @param boundary_value Value for Dirichlet or gradient for Neumann BC
 * @param boundary_face Face identifier (0=left, 1=right, 2=bottom, 3=top, 4=front, 5=back)
 * @param nx Number of grid points in x-direction
 * @param ny Number of grid points in y-direction
 * @param nz Number of grid points in z-direction
 * 
 * Applies the specified boundary condition to the given field.
 */
__global__ void applyBoundaryConditions(double* field_d,
                                       int boundary_type, double boundary_value,
                                       int boundary_face,
                                       int nx, int ny, int nz);

/**
 * @brief CUDA kernel for Jacobi iteration (for pressure Poisson equation)
 * @param p_d Device pointer to pressure field (input/output)
 * @param p_new_d Device pointer to new pressure field
 * @param rhs_d Device pointer to right-hand side
 * @param nx Number of grid points in x-direction
 * @param ny Number of grid points in y-direction
 * @param nz Number of grid points in z-direction
 * @param dx Grid spacing in x-direction
 * @param dy Grid spacing in y-direction
 * @param dz Grid spacing in z-direction
 * 
 * Performs one Jacobi iteration for solving the pressure Poisson equation.
 */
__global__ void jacobiIteration(const double* p_d, double* p_new_d, const double* rhs_d,
                               int nx, int ny, int nz,
                               double dx, double dy, double dz);

/**
 * @brief CUDA kernel for computing residuals
 * @param field_d Device pointer to field
 * @param field_old_d Device pointer to previous field values
 * @param residual_d Device pointer to output residual array
 * @param n Total number of grid points
 * 
 * Computes the L2 norm of the residual between current and previous
 * field values.
 */
__global__ void computeResidual(const double* field_d, const double* field_old_d,
                               double* residual_d, int n);

/**
 * @brief CUDA kernel for time integration (Runge-Kutta)
 * @param field_d Device pointer to field (input/output)
 * @param k1_d Device pointer to k1 coefficient
 * @param k2_d Device pointer to k2 coefficient
 * @param k3_d Device pointer to k3 coefficient
 * @param k4_d Device pointer to k4 coefficient
 * @param dt Time step size
 * @param n Total number of grid points
 * @param rk_stage Runge-Kutta stage (1-4)
 * 
 * Performs time integration using the 4th-order Runge-Kutta method.
 */
__global__ void rungeKuttaIntegration(double* field_d,
                                     const double* k1_d, const double* k2_d,
                                     const double* k3_d, const double* k4_d,
                                     double dt, int n, int rk_stage);

/**
 * @brief CUDA kernel for copying and scaling arrays
 * @param src_d Source array
 * @param dst_d Destination array
 * @param scale Scaling factor
 * @param n Array size
 * 
 * Utility kernel for copying arrays with optional scaling:
 * dst[i] = scale * src[i]
 */
__global__ void copyAndScale(const double* src_d, double* dst_d, double scale, int n);

/**
 * @brief CUDA kernel for initializing arrays with constant values
 * @param array_d Device pointer to array
 * @param value Initialization value
 * @param n Array size
 * 
 * Sets all elements of the array to the specified value.
 */
__global__ void initializeArray(double* array_d, double value, int n);

/**
 * @brief Host function to launch divergence computation kernel
 * @param u_d Device pointer to x-velocity
 * @param v_d Device pointer to y-velocity
 * @param w_d Device pointer to z-velocity
 * @param div_d Device pointer to divergence output
 * @param nx,ny,nz Grid dimensions
 * @param dx,dy,dz Grid spacing
 * @param stream CUDA stream for asynchronous execution
 * @return cudaError_t Error code
 */
cudaError_t launchDivergenceKernel(const double* u_d, const double* v_d, const double* w_d,
                                  double* div_d,
                                  int nx, int ny, int nz,
                                  double dx, double dy, double dz,
                                  cudaStream_t stream = 0);

/**
 * @brief Host function to launch pressure gradient kernel
 * @param p_d Device pointer to pressure
 * @param grad_px_d,grad_py_d,grad_pz_d Device pointers to gradient components
 * @param nx,ny,nz Grid dimensions
 * @param dx,dy,dz Grid spacing
 * @param stream CUDA stream for asynchronous execution
 * @return cudaError_t Error code
 */
cudaError_t launchPressureGradientKernel(const double* p_d,
                                        double* grad_px_d, double* grad_py_d, double* grad_pz_d,
                                        int nx, int ny, int nz,
                                        double dx, double dy, double dz,
                                        cudaStream_t stream = 0);

} // namespace kernels
} // namespace cfd

#endif // CUDA_KERNELS_CUH
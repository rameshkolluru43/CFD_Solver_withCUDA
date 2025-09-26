/**
 * @file cfd_solver.h
 * @brief Main CFD solver class definition with CUDA acceleration
 * @author Ramesh Kolluru
 * @date 2024
 * @version 0.1.0-alpha
 * 
 * This file contains the main CFD solver class that provides the interface
 * for computational fluid dynamics simulations using CUDA kernels for
 * GPU acceleration.
 */

#ifndef CFD_SOLVER_H
#define CFD_SOLVER_H

#include <vector>
#include <memory>
#include <string>
#include <cuda_runtime.h>

/**
 * @namespace cfd
 * @brief Main namespace for CFD solver components
 */
namespace cfd {

/**
 * @enum SolverType
 * @brief Enumeration of available solver types
 */
enum class SolverType {
    FINITE_VOLUME,    ///< Finite volume method
    FINITE_ELEMENT,   ///< Finite element method
    FINITE_DIFFERENCE ///< Finite difference method
};

/**
 * @enum TimeScheme
 * @brief Time integration schemes
 */
enum class TimeScheme {
    EXPLICIT_EULER,   ///< First-order explicit Euler
    RUNGE_KUTTA_2,    ///< Second-order Runge-Kutta
    RUNGE_KUTTA_4,    ///< Fourth-order Runge-Kutta
    IMPLICIT_EULER    ///< First-order implicit Euler
};

/**
 * @struct FluidProperties
 * @brief Structure containing fluid physical properties
 */
struct FluidProperties {
    double density;        ///< Fluid density [kg/m³]
    double viscosity;      ///< Dynamic viscosity [Pa·s]
    double specific_heat;  ///< Specific heat capacity [J/(kg·K)]
    double thermal_conductivity; ///< Thermal conductivity [W/(m·K)]
    
    /**
     * @brief Default constructor with air properties at STP
     */
    FluidProperties() : 
        density(1.225), 
        viscosity(1.81e-5), 
        specific_heat(1005.0),
        thermal_conductivity(0.0262) {}
};

/**
 * @struct SimulationParameters
 * @brief Configuration parameters for the simulation
 */
struct SimulationParameters {
    int nx;                ///< Number of grid points in x-direction
    int ny;                ///< Number of grid points in y-direction
    int nz;                ///< Number of grid points in z-direction
    double dx;             ///< Grid spacing in x-direction [m]
    double dy;             ///< Grid spacing in y-direction [m]
    double dz;             ///< Grid spacing in z-direction [m]
    double dt;             ///< Time step size [s]
    int max_iterations;    ///< Maximum number of time steps
    double convergence_tolerance; ///< Convergence criterion
    SolverType solver_type;    ///< Numerical solver type
    TimeScheme time_scheme;    ///< Time integration scheme
    FluidProperties fluid;     ///< Fluid properties
    
    /**
     * @brief Default constructor with typical values
     */
    SimulationParameters() :
        nx(128), ny(128), nz(1),
        dx(0.01), dy(0.01), dz(0.01),
        dt(1e-4), max_iterations(10000),
        convergence_tolerance(1e-6),
        solver_type(SolverType::FINITE_VOLUME),
        time_scheme(TimeScheme::RUNGE_KUTTA_4) {}
};

/**
 * @class CFDSolver
 * @brief Main CFD solver class with CUDA acceleration
 * 
 * This class provides the main interface for running computational fluid
 * dynamics simulations. It manages memory allocation on both CPU and GPU,
 * orchestrates the numerical solution process, and provides methods for
 * initialization, execution, and result extraction.
 * 
 * @example
 * ```cpp
 * cfd::CFDSolver solver;
 * cfd::SimulationParameters params;
 * params.nx = 256;
 * params.ny = 256;
 * solver.initialize(params);
 * solver.run();
 * auto results = solver.getResults();
 * ```
 */
class CFDSolver {
public:
    /**
     * @brief Default constructor
     */
    CFDSolver();
    
    /**
     * @brief Destructor - ensures proper cleanup of GPU memory
     */
    ~CFDSolver();
    
    /**
     * @brief Copy constructor (deleted to prevent accidental copying)
     */
    CFDSolver(const CFDSolver&) = delete;
    
    /**
     * @brief Assignment operator (deleted to prevent accidental copying)
     */
    CFDSolver& operator=(const CFDSolver&) = delete;
    
    /**
     * @brief Initialize the solver with given parameters
     * @param params Simulation parameters
     * @return true if initialization successful, false otherwise
     * 
     * This method allocates memory on both CPU and GPU, initializes
     * the computational grid, and sets up boundary conditions.
     */
    bool initialize(const SimulationParameters& params);
    
    /**
     * @brief Set boundary conditions
     * @param boundary_type Type of boundary condition
     * @param values Boundary values
     * @return true if successful, false otherwise
     */
    bool setBoundaryConditions(const std::string& boundary_type, 
                              const std::vector<double>& values);
    
    /**
     * @brief Set initial conditions
     * @param field_name Name of the field (e.g., "velocity", "pressure")
     * @param values Initial field values
     * @return true if successful, false otherwise
     */
    bool setInitialConditions(const std::string& field_name,
                             const std::vector<double>& values);
    
    /**
     * @brief Run the simulation
     * @return true if simulation completed successfully, false otherwise
     * 
     * This method executes the main time-stepping loop, calling the
     * appropriate CUDA kernels for each time step.
     */
    bool run();
    
    /**
     * @brief Run a single time step
     * @return true if time step completed successfully, false otherwise
     */
    bool timeStep();
    
    /**
     * @brief Check if simulation has converged
     * @return true if converged, false otherwise
     */
    bool hasConverged() const;
    
    /**
     * @brief Get current simulation time
     * @return Current time [s]
     */
    double getCurrentTime() const;
    
    /**
     * @brief Get current iteration number
     * @return Current iteration
     */
    int getCurrentIteration() const;
    
    /**
     * @brief Get residual values for monitoring convergence
     * @return Vector of residual values for each equation
     */
    std::vector<double> getResiduals() const;
    
    /**
     * @brief Extract simulation results
     * @param field_name Name of the field to extract
     * @return Vector containing field values
     */
    std::vector<double> getResults(const std::string& field_name = "pressure") const;
    
    /**
     * @brief Save results to file
     * @param filename Output filename
     * @param format File format ("vtk", "hdf5", "csv")
     * @return true if save successful, false otherwise
     */
    bool saveResults(const std::string& filename, const std::string& format = "vtk") const;
    
    /**
     * @brief Get memory usage information
     * @return Memory usage in bytes
     */
    size_t getMemoryUsage() const;
    
    /**
     * @brief Get CUDA device information
     * @return String containing device information
     */
    std::string getDeviceInfo() const;

private:
    struct Impl; ///< Forward declaration for PIMPL pattern
    std::unique_ptr<Impl> pImpl; ///< Private implementation pointer
};

/**
 * @brief Utility function to check CUDA errors
 * @param result CUDA runtime result
 * @param function Function name where error occurred
 * @param file Source file name
 * @param line Line number
 * @return true if no error, false otherwise
 */
bool checkCudaError(cudaError_t result, const char* function, 
                   const char* file, int line);

/**
 * @brief Macro for CUDA error checking
 */
#define CUDA_CHECK(call) checkCudaError((call), #call, __FILE__, __LINE__)

} // namespace cfd

#endif // CFD_SOLVER_H
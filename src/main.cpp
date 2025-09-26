/**
 * @file main.cpp
 * @brief Main entry point for the CFD solver application
 * @author Ramesh Kolluru
 * @date 2024
 * @version 0.1.0-alpha
 * 
 * This file contains the main function that serves as the entry point
 * for the CFD solver application. It handles command-line arguments,
 * initializes the solver, and manages the simulation execution.
 */

#include <iostream>
#include <string>
#include <chrono>
#include <memory>
#include <getopt.h>

#include "cfd_solver.h"

/**
 * @brief Print usage information
 * @param program_name Name of the executable
 */
void printUsage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [OPTIONS]\n"
              << "\nOptions:\n"
              << "  -h, --help              Show this help message\n"
              << "  -c, --config FILE       Configuration file path\n"
              << "  -i, --input FILE        Input mesh file\n"
              << "  -o, --output DIR        Output directory\n"
              << "  -t, --timesteps NUM     Number of time steps\n"
              << "  -d, --device ID         CUDA device ID to use\n"
              << "  -v, --verbose           Enable verbose output\n"
              << "  --nx NUM                Grid points in x-direction\n"
              << "  --ny NUM                Grid points in y-direction\n"
              << "  --nz NUM                Grid points in z-direction\n"
              << "  --dx FLOAT              Grid spacing in x-direction\n"
              << "  --dy FLOAT              Grid spacing in y-direction\n"
              << "  --dz FLOAT              Grid spacing in z-direction\n"
              << "  --dt FLOAT              Time step size\n"
              << "  --reynolds FLOAT        Reynolds number\n"
              << "\nExamples:\n"
              << "  " << program_name << " --config cavity_flow.cfg\n"
              << "  " << program_name << " --nx 256 --ny 256 --timesteps 1000\n"
              << "  " << program_name << " --input mesh.dat --output results/\n";
}

/**
 * @brief Parse command line arguments
 * @param argc Argument count
 * @param argv Argument vector
 * @param params Reference to simulation parameters to be filled
 * @param config_file Reference to configuration file path
 * @param input_file Reference to input file path
 * @param output_dir Reference to output directory path
 * @param device_id Reference to CUDA device ID
 * @param verbose Reference to verbose flag
 * @return true if parsing successful, false otherwise
 */
bool parseCommandLine(int argc, char* argv[], 
                     cfd::SimulationParameters& params,
                     std::string& config_file,
                     std::string& input_file,
                     std::string& output_dir,
                     int& device_id,
                     bool& verbose) {
    
    const struct option long_options[] = {
        {"help",      no_argument,       0, 'h'},
        {"config",    required_argument, 0, 'c'},
        {"input",     required_argument, 0, 'i'},
        {"output",    required_argument, 0, 'o'},
        {"timesteps", required_argument, 0, 't'},
        {"device",    required_argument, 0, 'd'},
        {"verbose",   no_argument,       0, 'v'},
        {"nx",        required_argument, 0, 1000},
        {"ny",        required_argument, 0, 1001},
        {"nz",        required_argument, 0, 1002},
        {"dx",        required_argument, 0, 1003},
        {"dy",        required_argument, 0, 1004},
        {"dz",        required_argument, 0, 1005},
        {"dt",        required_argument, 0, 1006},
        {"reynolds",  required_argument, 0, 1007},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "hc:i:o:t:d:v", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                printUsage(argv[0]);
                return false;
            case 'c':
                config_file = optarg;
                break;
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 't':
                params.max_iterations = std::stoi(optarg);
                break;
            case 'd':
                device_id = std::stoi(optarg);
                break;
            case 'v':
                verbose = true;
                break;
            case 1000:  // --nx
                params.nx = std::stoi(optarg);
                break;
            case 1001:  // --ny
                params.ny = std::stoi(optarg);
                break;
            case 1002:  // --nz
                params.nz = std::stoi(optarg);
                break;
            case 1003:  // --dx
                params.dx = std::stod(optarg);
                break;
            case 1004:  // --dy
                params.dy = std::stod(optarg);
                break;
            case 1005:  // --dz
                params.dz = std::stod(optarg);
                break;
            case 1006:  // --dt
                params.dt = std::stod(optarg);
                break;
            case 1007:  // --reynolds
                // Calculate viscosity from Reynolds number
                // Re = rho * U * L / mu
                // Assuming characteristic length L = dx and velocity U = 1.0
                params.fluid.viscosity = params.fluid.density * 1.0 * params.dx / std::stod(optarg);
                break;
            case '?':
                std::cerr << "Unknown option or missing argument\n";
                return false;
            default:
                return false;
        }
    }
    
    return true;
}

/**
 * @brief Load configuration from file
 * @param filename Configuration file path
 * @param params Reference to simulation parameters
 * @return true if successful, false otherwise
 */
bool loadConfiguration(const std::string& filename, cfd::SimulationParameters& params) {
    // TODO: Implement configuration file parsing
    // This is a placeholder for future implementation
    std::cout << "Loading configuration from: " << filename << std::endl;
    return true;
}

/**
 * @brief Main function
 * @param argc Argument count
 * @param argv Argument vector
 * @return Exit status
 */
int main(int argc, char* argv[]) {
    try {
        // Print banner
        std::cout << "CFD Solver with CUDA v0.1.0-alpha\n"
                  << "Copyright (C) 2024 Ramesh Kolluru\n"
                  << "This is free software; see the source for copying conditions.\n\n";
        
        // Initialize parameters with defaults
        cfd::SimulationParameters params;
        std::string config_file;
        std::string input_file;
        std::string output_dir = "output/";
        int device_id = 0;
        bool verbose = false;
        
        // Parse command line arguments
        if (!parseCommandLine(argc, argv, params, config_file, input_file, 
                             output_dir, device_id, verbose)) {
            return 1;
        }
        
        // Load configuration file if specified
        if (!config_file.empty()) {
            if (!loadConfiguration(config_file, params)) {
                std::cerr << "Error: Failed to load configuration file: " << config_file << std::endl;
                return 1;
            }
        }
        
        // Set CUDA device
        cudaError_t cuda_status = cudaSetDevice(device_id);
        if (cuda_status != cudaSuccess) {
            std::cerr << "Error: Failed to set CUDA device " << device_id 
                      << ": " << cudaGetErrorString(cuda_status) << std::endl;
            return 1;
        }
        
        if (verbose) {
            std::cout << "Using CUDA device: " << device_id << std::endl;
        }
        
        // Create and initialize solver
        std::unique_ptr<cfd::CFDSolver> solver = std::make_unique<cfd::CFDSolver>();
        
        if (verbose) {
            std::cout << solver->getDeviceInfo() << std::endl;
            std::cout << "Grid dimensions: " << params.nx << " x " << params.ny << " x " << params.nz << std::endl;
            std::cout << "Grid spacing: " << params.dx << " x " << params.dy << " x " << params.dz << std::endl;
            std::cout << "Time step: " << params.dt << std::endl;
            std::cout << "Max iterations: " << params.max_iterations << std::endl;
        }
        
        // Initialize solver
        auto start_init = std::chrono::high_resolution_clock::now();
        if (!solver->initialize(params)) {
            std::cerr << "Error: Failed to initialize solver" << std::endl;
            return 1;
        }
        auto end_init = std::chrono::high_resolution_clock::now();
        
        if (verbose) {
            auto init_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_init - start_init);
            std::cout << "Initialization completed in " << init_time.count() << " ms" << std::endl;
            std::cout << "Memory usage: " << solver->getMemoryUsage() / (1024*1024) << " MB" << std::endl;
        }
        
        // Run simulation
        std::cout << "Starting simulation..." << std::endl;
        auto start_sim = std::chrono::high_resolution_clock::now();
        
        bool success = solver->run();
        
        auto end_sim = std::chrono::high_resolution_clock::now();
        auto sim_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_sim - start_sim);
        
        if (success) {
            std::cout << "Simulation completed successfully in " << sim_time.count() << " ms" << std::endl;
            
            // Save results
            std::string result_file = output_dir + "/results.vtk";
            if (solver->saveResults(result_file)) {
                std::cout << "Results saved to: " << result_file << std::endl;
            } else {
                std::cerr << "Warning: Failed to save results" << std::endl;
            }
            
            if (verbose) {
                auto residuals = solver->getResiduals();
                std::cout << "Final residuals: ";
                for (size_t i = 0; i < residuals.size(); ++i) {
                    std::cout << residuals[i];
                    if (i < residuals.size() - 1) std::cout << ", ";
                }
                std::cout << std::endl;
            }
        } else {
            std::cerr << "Error: Simulation failed" << std::endl;
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Error: Unknown exception occurred" << std::endl;
        return 1;
    }
    
    return 0;
}
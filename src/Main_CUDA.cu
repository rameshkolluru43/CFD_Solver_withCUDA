#include <cuda_runtime.h>
#include <iostream>
#include <string>
#include <vector>

extern std::string Test_Case_JSON_File;
extern std::string Error_File;
extern std::string Solution_File;
extern int Test_Case;
extern bool Is_Viscous_Wall;

void readJSON(const std::string &);
void parseTestCaseJSON(const std::string &);
bool initializeGridFiles();
void createOutputDirectories();
void Initialize_TestCase();
bool Inviscid_Solver(std::string &, std::string &);
bool Viscous_Solver(std::string &, std::string &);

static void print_gpu_info()
{
    int deviceCount = 0;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        std::cout << "[CUDA] No CUDA-capable devices found." << std::endl;
        return;
    }
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "\n=== CUDA Device ===" << std::endl;
    std::cout << "  Device: " << prop.name << std::endl;
    std::cout << "  Compute Capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "  Global Memory: " << prop.totalGlobalMem / (1024*1024) << " MB" << std::endl;
    std::cout << "  SM Count: " << prop.multiProcessorCount << std::endl;
    std::cout << "===================" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.json>" << std::endl;
        return EXIT_FAILURE;
    }

    print_gpu_info();

    std::string jsonFile = argv[1];
    try {
        readJSON(jsonFile);
    } catch (const std::exception &e) {
        std::cerr << "Error reading JSON: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    if (Test_Case_JSON_File.empty()) {
        std::cerr << "Error: Test_Case_JSON_File not set." << std::endl;
        return EXIT_FAILURE;
    }

    parseTestCaseJSON(Test_Case_JSON_File);

    if (!initializeGridFiles()) {
        std::cerr << "Error: Failed to initialize grid files." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Checking and creating the directories for solution files" << std::endl;
    createOutputDirectories();
    Initialize_TestCase();

    std::cout << "*******----------GPU-Enabled Solver Running-----------**********\n";

    bool solver_success = false;
    if (Is_Viscous_Wall) {
        std::cout << "Viscous Terms are Enabled solving NS Solver" << std::endl;
        solver_success = Viscous_Solver(Error_File, Solution_File);
    } else {
        std::cout << "Inviscid Solver ---- Solving Euler Equations" << std::endl;
        solver_success = Inviscid_Solver(Error_File, Solution_File);
    }

    if (!solver_success) {
        std::cerr << "Error: Solver failed." << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "GPU-enabled CFD simulation completed successfully!" << std::endl;
    return EXIT_SUCCESS;
}

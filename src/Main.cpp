/**
 * @file Main.cpp
 * @brief Main entry point for the simulation program.
 * Contains the main function and initialization routines.
 * Reads the input JSON file and initializes the simulation parameters.
 * Calls the appropriate solver based on the configuration.
 */
#include "definitions.h"
#include "Globals.h"
#include "Main.h"
#include "Directory_Files.h"
#include "Solver.h"
#include "MPI_Utils.h"

void readInputFile(int argc, char *argv[])
{
    if (argc < 2)
    {
        cerr << "Error: No JSON input file provided." << endl;
        exit(EXIT_FAILURE);
    }
    string jsonFileName = argv[1];
    try
    {
        readJSON(jsonFileName);
    }
    catch (const std::exception &e)
    {
        cerr << "Error while reading JSON file: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    if (Test_Case_JSON_File.empty())
    {
        cerr << "Error: Test_Case_JSON_File is not set. Please check the JSON file or initialization." << endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Executes the appropriate solver based on the simulation configuration.
 *
 * This function validates the test case setting and then dispatches control to the corresponding
 * solver:
 * - If the test case is not set (i.e., Test_Case equals 0), it logs an error message and returns false.
 * - If the viscous wall condition (Is_Viscous_Wall) is active, it runs the Navier-Stokes solver with viscous terms.
 * - Otherwise, it runs the solver for Euler equations.
 *
 * @details
 * The function verifies that Test_Case is properly set. A Test_Case value of 0 is considered an error since test case numbers
 * should start at 1. Upon validation, it prints a message indicating the current solver in use, then calls either the
 * Viscous_Solver or Inviscid_Solver function based on the Is_Viscous_Wall flag.
 *
 * @note Proper initialization is assumed to be performed elsewhere; the Test_Case should be assigned a valid, non-zero value.
 *
 * @return true if solver completes successfully, false otherwise
 *
 * @see Viscous_Solver, Inviscid_Solver
 */
bool runSolver()
{
    try
    {
        if (Test_Case == 0)
        {
            cerr << "Error: Test_Case is not set. Please check the JSON file or initialization." << endl;
            cerr << "Test Case Number starts from 1\n";
            return false;
        }

        if (CFD_MPI_Is_Root())
            cout << "*******----------The Following Solver is currently Running-----------**********\n";

        bool solver_success = false;
        if (Is_Viscous_Wall)
        {
            if (CFD_MPI_Is_Root())
                cout << "Viscous Terms are Enabled solving NS Solver" << endl;
            solver_success = Viscous_Solver(Error_File, Solution_File);
            if (!solver_success)
            {
                cerr << "Error: Viscous solver failed to complete successfully" << endl;
                return false;
            }
        }
        else
        {
            if (CFD_MPI_Is_Root())
                cout << "Inviscid Solver ---- Solving Euler Equations" << endl;
            solver_success = Inviscid_Solver(Error_File, Solution_File);
            if (!solver_success)
            {
                cerr << "Error: Inviscid solver failed to complete successfully" << endl;
                return false;
            }
        }

        if (CFD_MPI_Is_Root())
            cout << "Solver completed successfully" << endl;
        return true;
    }
    catch (const std::exception &e)
    {
        cerr << "Exception in runSolver: " << e.what() << endl;
        return false;
    }
    catch (...)
    {
        cerr << "Unknown exception occurred in runSolver" << endl;
        return false;
    }
}

int main(int argc, char *argv[])
{
    CFD_MPI_Initialize(&argc, &argv);
    if (CFD_MPI_Is_Enabled() && CFD_MPI_Is_Root())
    {
        cout << "MPI enabled with " << CFD_MPI_Size() << " ranks" << endl;
    }

    // Read input JSON file and initialize parameters
    readInputFile(argc, argv);
    // Reads the test case JSON file and populates the parameters
    parseTestCaseJSON(Test_Case_JSON_File);

    if (!initializeGridFiles())
    {
        cerr << "Error: Failed to initialize grid files. Exiting..." << endl;
        CFD_MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    if (CFD_MPI_Is_Root())
    {
        cout << "Checking and creating the directories for solution files" << endl;
        createOutputDirectories();
    }
    CFD_MPI_Barrier();
    Initialize_TestCase();

    if (!runSolver())
    {
        cerr << "Error: Solver execution failed. Exiting..." << endl;
        CFD_MPI_Finalize();
        return EXIT_FAILURE;
    }

    if (CFD_MPI_Is_Root())
        cout << "CFD simulation completed successfully!" << endl;
    CFD_MPI_Finalize();
    return EXIT_SUCCESS;
}

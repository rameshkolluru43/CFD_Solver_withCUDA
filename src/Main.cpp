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

void initializeGridDimensions()
{
    nx_1 = 0;
    ny_1 = 0;
}

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
 * - If the test case is not set (i.e., Test_Case equals 0), it logs an error message and terminates the program.
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
 * @see Viscous_Solver, Inviscid_Solver
 */
void runSolver()
{
    if (Test_Case == 0)
    {
        cerr << "Error: Test_Case is not set. Please check the JSON file or initialization." << endl;
        cerr << "Test Case Number starts from 1\n";
        exit(EXIT_FAILURE);
    }
    cout << "*******----------The Following Solver is currently Running-----------**********\n";
    if (Is_Viscous_Wall)
    {
        cout << "Viscous Terms are Enabled solving NS Solver" << endl;
        Viscous_Solver(Error_File, Solution_File);
    }
    else
    {
        cout << "Inviscid Solver ---- Solving Euler Equations" << endl;
        Inviscid_Solver(Error_File, Solution_File);
    }
}

int main(int argc, char *argv[])
{
    initializeGridDimensions();
    readInputFile(argc, argv);
    readTestCaseJSON(Test_Case_JSON_File);
    runSolver();
    return 0;
}

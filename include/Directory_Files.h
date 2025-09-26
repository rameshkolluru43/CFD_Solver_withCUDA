#ifndef DIRECTORY_FILES_H
#define DIRECTORY_FILES_H
#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"
#include "Boundary_Conditions.h"

//------------Functions for creating Test case directories and reading JSON files----------------
void CreateTestCaseDirectories(int &);
void readJSON(const string &);
void testCase(int &);
void searchGridFiles();
void findGridFiles(string &, int &, int &, vector<string> &);
void Directory_Name();
void File_Name();
void readTestCaseJSON(const std::string &);
bool readBoundaryConditions(const std::string &, InletCondition &, ExitCondition &);
bool readInitialConditionsjson(const std::string &, InitialCondition &);
bool readInitialAndBoundaryConditions(const std::string &, InitialCondition &, InletCondition &, ExitCondition &);
void initializeGridFiles();
#endif // #ifndef DIRECTORY_FILES_H
#ifndef DIRECTORY_FILES_H
#define DIRECTORY_FILES_H
#include "definitions.h"
#include "Globals.h"
#include "Test_Cases.h"
#include "Boundary_Conditions.h"

//------------Functions for creating Test case directories and reading JSON files----------------
void CreateTestCaseDirectories(int &);
void createOutputDirectories();
void readJSON(const string &);
void searchGridFiles();
void findGridFiles(string &, int &, int &, vector<string> &);
void Directory_Name();
void File_Name();
void parseTestCaseJSON(const std::string &);
bool readBoundaryConditions(const std::string &, InletCondition &, ExitCondition &);
bool readInitialConditionsjson(const std::string &, InitialCondition &);
bool readInitialAndBoundaryConditions(const std::string &, InitialCondition &, InletCondition &, ExitCondition &);
bool initializeGridFiles();
#endif // #ifndef DIRECTORY_FILES_H
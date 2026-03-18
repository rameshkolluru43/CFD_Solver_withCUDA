#include "definitions.h"
#include "Globals.h"
#include "IO_Write.h"

void createOutputDirectories()
{
    // Base directory for results
    string baseDir = "../2D_Euler_Solutions/";
    // Check if the directory already exists, if not, create it
    if (!filesystem::exists(baseDir))
    {
        cout << "Base Directory: " << baseDir << endl;
        filesystem::create_directories(baseDir);
    }
    else
    {
        cerr << "Error: Directory already exists!" << endl;
        //        return;
    }

    // Create directory based on test case
    string testCaseDir = baseDir + Test_Case_Name;
    if (!filesystem::exists(testCaseDir))
    {
        cout << "Test Case Directory: " << testCaseDir << endl;
        filesystem::create_directories(testCaseDir);
    }
    else
    {
        cerr << "Error: Directory already exists!" << endl;
        //      return;
    }

    // Create directories based on numerical method
    string methodDir = testCaseDir + "/" + (Is_Implicit_Method ? "Implicit" : "Explicit");
    if (!filesystem::exists(methodDir))
    {
        cout << "Method Directory: " << methodDir << endl;
        filesystem::create_directories(methodDir);
    }
    else
    {
        cerr << "Error: Directory already exists!" << endl;
        // return;
    }

    // Create flux-based directory
    string fluxDir = methodDir + "/Flux_" + to_string(Flux_Type);
    if (!filesystem::exists(fluxDir))
    {
        cout << "Flux Type Directory : " << fluxDir << endl;
        filesystem::create_directories(fluxDir);
    }
    else
    {
        cerr << "Error: Directory already exists!" << endl;
        // return;
    }

    // Create grid-specific directory
    string gridDir = fluxDir + "/GridSize_" + to_string(Grid_Size);
    if (!filesystem::exists(gridDir))
    {
        cout << "Grid Size Directory: " << gridDir << endl;
        filesystem::create_directories(gridDir);
    }
    else
    {
        cerr << "Error: Directory already exists!" << endl;
        //  return;
    }

    // Create additional directory for WENO if enabled
    if (Is_WENO)
    {
        string wenoDir = gridDir + "/WENO";
        if (!filesystem::exists(wenoDir))
        {
            filesystem::create_directories(wenoDir);
        }
        else
        {
            cerr << "Error: Directory already exists!" << endl;
            //     return;
        }
        gridDir = wenoDir;
    }

    // Create output files in the final directory

    string outputFilePath = gridDir + "/results.txt";
    Solution_File = gridDir + "/Solution" + to_string(Flux_Type) + to_string(Grid_Size) + ".txt";
    Error_File = gridDir + "/Error.txt";
    Initial_Solution_File = gridDir + "/Initial_Solution.txt";
    Final_Solution_File = gridDir + "/Final_Solution.vtk";

    cout << Solution_File << endl;
    cout << Error_File << endl;
    cout << Initial_Solution_File << endl;
    cout << Final_Solution_File << endl;

    ofstream outputFile(outputFilePath);
    if (outputFile.is_open())
    {
        outputFile << "Simulation Results for Test Case: " << Test_Case << "\n";
        outputFile << "Flux Type: " << Flux_Type << "\n";
        outputFile << "Numerical Method: " << (Is_Implicit_Method ? "Implicit" : "Explicit") << "\n";
        outputFile << "Grid Size: " << Grid_Size << "\n";
        outputFile << "WENO Enabled: " << (Is_WENO ? "Yes" : "No") << "\n";
        outputFile.close();
        cout << "Output written to: " << outputFilePath << endl;
    }
    else
    {
        cerr << "Error: Unable to create output file!" << endl;
    }
}

void findGridFiles(string &folderPath, int &desiredX, int &desiredY, vector<string> &gridFiles)
{

    // string pattern = Test_Case_Name + "_" + to_string(desiredX) + "_" + to_string(desiredY) + ".vtk";
    string pattern = Test_Case_Name + R"((\d+)_(\d+).vtk)";
    boost::regex filePattern(pattern);
    boost::smatch match;
    cout << "Folder Path: " << folderPath << endl;
    cout << "Pattern: " << pattern << endl;
    // findGridFiles Function :
    // • The function takes a folder path and desired grid dimensions.
    // • It uses a regular expression(Flow_Over_Cylinder_(\d +) _(\d +)\.txt) to match and capture the grid sizes from each filename.
    // • It compares the captured grid sizes to the desired values and, if they match, adds the file’s full path to the results.return gridFiles;
    for (const auto &entry : filesystem::directory_iterator(folderPath))
    {
        string fileName = entry.path().string();
        if (boost::regex_search(fileName, match, filePattern))
        {
            int x = std::stoi(match[1].str());
            int y = std::stoi(match[2].str());
            cout << "X: " << x << " Y: " << y << endl;
            if (x == desiredX && y == desiredY)
            {
                gridFiles.push_back(entry.path().string());
            }
            else
            {
                cout << "Could not find grid file for desired dimensions ..... exiting" << endl;
                exit(0);
            }
        }
    }
}

void searchGridFiles()
{
    string gridFile = "../Grid_Files/TestCase_" + to_string(Test_Case) + "/GridSize_" + to_string(Grid_Size) + ".vtk";
    cout << "Searching for grid file: " << gridFile << endl;
    if (filesystem::exists(gridFile))
    {
        Grid_File = gridFile;
    }
    else
    {
        cerr << "Error: Grid file not found!" << endl;
    }
}
#include "definitions.h"
#include "Globals.h"
#include "IO_Write.h"

static string getLimiterName(int limiterCase);

static string getDissipationName(int type)
{
    switch (type)
    {
    case 1:
        return "LLF";
    case 2:
        return "MOVERS";
    case 3:
        return "ROE";
    case 4:
        return "RICCA";
    case 5:
        return "MOVERS_NWSC";
    case 6:
        return "RICCA_LLF";
    default:
        return "UNKNOWN_" + to_string(type);
    }
}

static string getFluxSchemeName()
{
    if (Is_WENO)
        return "WENO";

    const string orderName = Is_Second_Order ? "2O" : "1O";
    return orderName + "_" + getLimiterName(Limiter_Case);
}

static string getLimiterName(int limiterCase)
{
    switch (limiterCase)
    {
    case 0:
        return "MinMod";
    case 1:
        return "Superbee";
    case 2:
        return "MCD";
    case 3:
        return "VanLeer";
    case 4:
        return "VanAlbada";
    case 5:
        return "Venkat";
    default:
        return "Limiter" + to_string(limiterCase);
    }
}

void createOutputDirectories()
{
    try
    {
        // Base directory for results
        string baseDir = "../solutions/2D_Euler_Solutions/";
        // Check if the directory already exists, if not, create it
        if (!filesystem::exists(baseDir))
        {
            cout << "Base Directory: " << baseDir << endl;
            filesystem::create_directories(baseDir);
        }
        else if (!filesystem::is_directory(baseDir))
        {
            throw runtime_error("Output base path exists but is not a directory: " + baseDir);
        }

        // Create directory based on test case
        string testCaseDir = baseDir + Test_Case_Name;
        if (!filesystem::exists(testCaseDir))
        {
            cout << "Test Case Directory: " << testCaseDir << endl;
            filesystem::create_directories(testCaseDir);
        }
        else if (!filesystem::is_directory(testCaseDir))
        {
            throw runtime_error("Test case output path exists but is not a directory: " + testCaseDir);
        }

        // Create directories based on numerical method
        string methodDir = testCaseDir + "/" + (Is_Implicit_Method ? "Implicit" : "Explicit");
        if (!filesystem::exists(methodDir))
        {
            cout << "Method Directory: " << methodDir << endl;
            filesystem::create_directories(methodDir);
        }
        else if (!filesystem::is_directory(methodDir))
        {
            throw runtime_error("Method output path exists but is not a directory: " + methodDir);
        }

        // Create flux-based directory
        string fluxDir = methodDir + "/Flux_" + to_string(Flux_Type);
        if (!filesystem::exists(fluxDir))
        {
            cout << "Flux Type Directory : " << fluxDir << endl;
            filesystem::create_directories(fluxDir);
        }
        else if (!filesystem::is_directory(fluxDir))
        {
            throw runtime_error("Flux output path exists but is not a directory: " + fluxDir);
        }

        // Create grid-specific directory
        string gridDir = fluxDir + "/GridSize_" + to_string(Grid_Size);
        if (!filesystem::exists(gridDir))
        {
            cout << "Grid Size Directory: " << gridDir << endl;
            filesystem::create_directories(gridDir);
        }
        else if (!filesystem::is_directory(gridDir))
        {
            throw runtime_error("Grid output path exists but is not a directory: " + gridDir);
        }

        // Create additional directory for WENO if enabled
        if (Is_WENO)
        {
            string wenoDir = gridDir + "/WENO";
            if (!filesystem::exists(wenoDir))
            {
                filesystem::create_directories(wenoDir);
            }
            else if (!filesystem::is_directory(wenoDir))
            {
                throw runtime_error("WENO output path exists but is not a directory: " + wenoDir);
            }
            gridDir = wenoDir;
        }

        // Create output files in the final directory with descriptive names
        string runTag;
        if (Is_WENO)
        {
            runTag = "WENO_Grid_Size_" + to_string(Grid_Size);
        }
        else
        {
            const string orderName = Is_Second_Order ? "2O" : "1O";
            runTag = getDissipationName(Dissipation_Type) + "_" + orderName + "_" + getLimiterName(Limiter_Case);
        }

        string outputFilePath = gridDir + "/results_" + runTag + ".txt";
        Solution_File = gridDir + "/Solution_" + runTag + ".txt";
        Error_File = gridDir + "/Error_" + runTag + ".txt";
        Initial_Solution_File = gridDir + "/Initial_Solution_" + runTag + ".txt";
        Final_Solution_File = gridDir + "/" + Test_Case_Name +
                              "_Grid_Size_" + to_string(Grid_Size) +
                              "_FluxScheme_" + getFluxSchemeName() +
                              "_Dissipation_" + getDissipationName(Dissipation_Type) +
                              ".vtk";
        CF_File = gridDir + "/CF_" + runTag + ".txt";
        QW_File = gridDir + "/QW_" + runTag + ".txt";

        cout << Solution_File << endl;
        cout << Error_File << endl;
        cout << Initial_Solution_File << endl;
        cout << Final_Solution_File << endl;
        if (Is_Viscous_Wall)
        {
            cout << CF_File << endl;
            cout << QW_File << endl;
        }

        ofstream outputFile = CFD_OpenOutputFileOrThrow(outputFilePath, "createOutputDirectories");
        outputFile << "Simulation Results for Test Case: " << Test_Case << "\n";
        outputFile << "Output folder Flux_Type (label only): " << Flux_Type << "\n";
        outputFile << "Dissipation_Type (actual face flux / dissipation): " << Dissipation_Type
                   << "  (1=LLF 2=MOVERS 3=Roe 4=RICCA 5=MOVERS_NWSC)\n";
        outputFile << "Numerical Method: " << (Is_Implicit_Method ? "Implicit" : "Explicit") << "\n";
        outputFile << "Grid Size: " << Grid_Size << "\n";
        outputFile << "WENO Enabled: " << (Is_WENO ? "Yes" : "No") << "\n";
        outputFile.close();
        cout << "Output written to: " << outputFilePath << endl;
    }
    catch (const filesystem::filesystem_error &e)
    {
        throw runtime_error(string("createOutputDirectories: filesystem error: ") + e.what());
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
    if (!filesystem::exists(folderPath) || !filesystem::is_directory(folderPath))
        throw runtime_error("findGridFiles: grid folder not found: " + folderPath);
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
                throw runtime_error("findGridFiles: grid file dimensions do not match requested size");
            }
        }
    }
    if (gridFiles.empty())
        throw runtime_error("findGridFiles: no matching grid files found in " + folderPath);
}

void searchGridFiles()
{
    string gridFile = "../input/Grid_Files/TestCase_" + to_string(Test_Case) + "/GridSize_" + to_string(Grid_Size) + ".vtk";
    cout << "Searching for grid file: " << gridFile << endl;
    if (filesystem::exists(gridFile))
    {
        Grid_File = gridFile;
    }
    else
    {
        throw runtime_error("searchGridFiles: grid file not found: " + gridFile);
    }
}
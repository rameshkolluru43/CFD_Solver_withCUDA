#include "definitions.h"
#include "Globals.h"
#include "Directory_Files.h"
#include "Test_Cases.h"

/// @brief
/// @param jsonFileName  The JSON file name
/// @return void     No return value
/// @details This function reads the JSON file and sets the parameters for the solver and all other parameters
double Limiter_Zeta, Limiter_Zeta1;
string Flow_Type, Test_Case_JSON_File, Test_Case_Config_File, Solver_Name, Description, Author, GeometryType;
bool Is_Conservative, Is_Viscous;
int Solver_Type;
InletCondition inletCond;
ExitCondition exitCond;
InitialCondition initCond;
WallCondition wallCond;
GeometryParams geomParams;
MeshParams meshParams;
vector<string> gridFiles;
string gridDir, Test_Case_Name = "Default_Test_Case", GridVTKFile;

void readJSON(const std::string &jsonFileName)
{
    std::ifstream fileStream(jsonFileName);
    if (!fileStream.is_open())
    {
        throw std::runtime_error("Error: JSON input file not found: " + jsonFileName);
    }

    Json::Reader reader;
    Json::Value root;
    if (!reader.parse(fileStream, root))
    {
        throw std::runtime_error("Error: Failed to parse JSON input.");
    }

    // Extract Solver Config First
    const auto &testCaseData = root["TestCase"];
    const auto &simData = root["Simulation"];
    const auto &solverData = root["Solver"];
    const auto &limiterData = root["LimiterCoefficients"];

    // Test case details
    Test_Case = testCaseData["Test_Case"].asInt(); // optional use
    Test_Case_Name = testCaseData["Test_Case_Name"].asString();
    Test_Case_JSON_File = testCaseData["Test_Case_Json"].asString();

    // Simulation settings
    Initialize_Type = simData["Initialize_Type"].asInt();
    Is_Implicit_Method = simData["Is_Implicit_Method"].asBool();
    Total_Iterations = simData["Total_Iterations"].asInt();
    CFL = simData["CFL"].asDouble();
    Terminating_Time = simData["Terminating_Time"].asDouble();
    Is_Time_Dependent = simData["Is_Time_Dependent"].asBool();

    // Solver settings
    Solver_Type = solverData["Solver_Type"].asInt();          // new
    Solver_Name = solverData["Solver_Name"].asString();       // new
    Is_Conservative = solverData["Is_Conservative"].asBool(); // new
    Is_Viscous = solverData["Is_Viscous"].asBool();           // new
    Limiter_Case = solverData["Limiter_Case"].asInt();
    Area_Weighted_Average = solverData["Area_Weighted_Average"].asInt();
    Flux_Type = solverData["Flux_Type"].asInt();
    NUM_FLUX_COMPONENTS = solverData["NUM_FLUX_COMPONENTS"].asInt();
    Is_Second_Order = solverData["Is_Second_Order"].asBool();
    Time_Accurate = solverData["Time_Accurate"].asBool();
    Local_Time_Stepping = solverData["Local_Time_Stepping"].asBool();
    Non_Dimensional_Form = solverData["Non_Dimensional_Form"].asBool();
    Is_WENO = solverData["Is_WENO"].asBool();
    Dissipation_Type = solverData["Dissipation_Type"].asInt();
    Is_MOVERS_1 = solverData["Is_MOVERS_1"].asBool();
    Enable_Entropy_Fix = solverData["Enable_Entropy_Fix"].asBool();
    if (solverData.isMember("Enable_AMR"))
        Enable_AMR = solverData["Enable_AMR"].asBool();
    if (solverData.isMember("AMR_Period"))
        AMR_Period = solverData["AMR_Period"].asInt();
    if (solverData.isMember("AMR_Gradient_Threshold"))
        AMR_Gradient_Threshold = solverData["AMR_Gradient_Threshold"].asDouble();
    if (solverData.isMember("AMR_Max_Fraction"))
        AMR_Max_Fraction = solverData["AMR_Max_Fraction"].asDouble();

    // Limiter coefficients
    Limiter_Zeta = limiterData["Limiter_Zeta"].asDouble();
    Limiter_Zeta1 = limiterData["Limiter_Zeta1"].asDouble();
    // Epsilon = limiterData["Epsilon"].asDouble(); // new
}

// Function to read parameters from a configuration file
map<string, string> ReadConfigFile(const string &filename)
{
    map<string, string> config;
    ifstream file(filename);

    if (!file.is_open())
    {
        cerr << "Error: Unable to open configuration file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    string line;
    while (getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue; // Skip empty lines and comments
        size_t delimiter = line.find('=');
        if (delimiter == string::npos)
            continue;

        string key = line.substr(0, delimiter);
        string value = line.substr(delimiter + 1);

        // Trim whitespace from key and value
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);

        config[key] = value;
    }
    file.close();
    return config;
}

// Reads the test case JSON file and initializes the grid and data structures
// from the grid file name initializes the grid and data structures and initializes the test case
void parseTestCaseJSON(const std::string &jsonFileName)
{
    std::ifstream fileStream(jsonFileName);
    if (!fileStream.is_open())
    {
        throw std::runtime_error("Error: Test Case JSON file not found: " + jsonFileName);
    }

    Json::Reader reader;
    Json::Value root;
    if (!reader.parse(fileStream, root))
    {
        throw std::runtime_error("Error: Failed to parse Test Case JSON.");
    }

    // Extract General Info
    const auto &genInfo = root["GeneralDetails"];
    Test_Case_Name = genInfo["name"].asString();
    // Description = genInfo["description"].asString();
    // GeometryType = genInfo["type"].asString();
    // Author = genInfo["author"].asString();

    // Geometry Parameters
    const auto &geom = root["geometry_parameters"];
    geomParams.radius = geom["radius"].asDouble();
    geomParams.length = geom["length"].asDouble();
    geomParams.thickness = geom["thickness"].asDouble();

    // Mesh Parameters
    const auto &mesh = root["mesh_parameters"];
    meshParams.gridSize = mesh["Grid_Size"].asInt();
    meshParams.nx = mesh["Nx"].asInt();
    meshParams.ny = mesh["Ny"].asInt();
    meshParams.meshType = mesh["mesh_type"].asString();

    // Flow Conditions
    const auto &flow = root["Flow_Conditions"];

    const auto &inlet = flow["inlet_conditions"];
    inletCond.type = inlet["inlet_type"].asString();
    inletCond.P = inlet["Pressure_Static_Inlet"].asDouble();
    inletCond.Rho = inlet["Rho_Static_Inlet"].asDouble();
    inletCond.M = inlet["Inlet_Mach_No"].asDouble();
    inletCond.u = inlet["u"].asDouble();
    inletCond.v = inlet["v"].asDouble();

    const auto &exit = flow["exit_conditions"];
    exitCond.type = exit["exit_type"].asString();
    exitCond.P = exit["Pressure_Static_Exit"].asDouble();
    exitCond.Rho = exit["Rho_Static_Exit"].asDouble();
    exitCond.M = exit["Exit_Mach_No"].asDouble();
    exitCond.u = exit["u"].asDouble();
    exitCond.v = exit["v"].asDouble();

    const auto &wall = flow["wall_conditions"];
    wallCond.type = wall["wall_type"].asString();
    wallCond.T = wall["wall_temperature"].asDouble();
    wallCond.u = wall["wall_velocity"].asDouble();
    wallCond.v = wall["wall_velocity"].asDouble();

    const auto &init = flow["Initial_Conditions"];
    initCond.P = init["Pressure_Static_Inlet"].asDouble();
    initCond.Rho = init["Rho_Static_Inlet"].asDouble();
    initCond.M = init["Inlet_Mach_No"].asDouble();
    initCond.u = init["u"].asDouble();
    initCond.v = init["v"].asDouble();
}

// Reads boundary conditions from a JSON file and populates the inlet and exit structures
bool readBoundaryConditionsjson(const std::string &filename, InletCondition &inletCond, ExitCondition &exitCond)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }
    Json::Reader reader;
    Json::Value BoundaryConditions;
    if (!reader.parse(file, BoundaryConditions))
    {
        throw std::runtime_error("Error: Failed to parse JSON input.");
    }

    // Extract inlet conditions
    inletCond.P = BoundaryConditions["inlet_conditions"]["Pressure_Static_Inlet"].asDouble();
    inletCond.Rho = BoundaryConditions["inlet_conditions"]["Rho_Static_Inlet"].asDouble();
    inletCond.M = BoundaryConditions["inlet_conditions"]["Inlet_Mach_No"].asDouble();
    inletCond.u = BoundaryConditions["inlet_conditions"]["V_1"].asDouble();
    inletCond.v = BoundaryConditions["inlet_conditions"]["V_2"].asDouble();
    inletCond.T = BoundaryConditions["inlet_conditions"]["Temperature_Static_Inlet"].asDouble();
    inletCond.T_inf = BoundaryConditions["inlet_conditions"]["T_inf"].asDouble();
    inletCond.P_inf = BoundaryConditions["inlet_conditions"]["P_inf"].asDouble();
    inletCond.Rho_inf = BoundaryConditions["inlet_conditions"]["Rho_inf"].asDouble();

    // Extract exit conditions
    exitCond.P = BoundaryConditions["exit_conditions"]["exit_pressure"].asDouble();
    exitCond.Rho = BoundaryConditions["exit_conditions"]["exit_density"].asDouble();
    exitCond.T = BoundaryConditions["exit_conditions"]["exit_temperature"].asDouble();
    exitCond.M = BoundaryConditions["exit_conditions"]["exit_mach_no"].asDouble();
    exitCond.u = BoundaryConditions["exit_conditions"]["V_1"].asDouble();
    exitCond.v = BoundaryConditions["exit_conditions"]["V_2"].asDouble();

    return true;
}

// Reads boundary conditions from a JSON file and populates the inlet and exit structures
bool readInitialConditionsjson(const std::string &filename, InitialCondition &initialCond)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }
    Json::Reader reader;
    Json::Value InitialConditions;
    if (!reader.parse(file, InitialConditions))
    {
        throw std::runtime_error("Error: Failed to parse JSON input.");
    }

    // Extract inlet conditions
    // initialCond.testCase = InitialConditions["initial_conditions"]["test_case"].asInt();
    initialCond.P = InitialConditions["Initial_Conditions"]["Pressure_Static_Inlet"].asDouble();
    initialCond.Rho = InitialConditions["Initial_Conditions"]["Rho_Static_Inlet"].asDouble();
    initialCond.M = InitialConditions["Initial_Conditions"]["Inlet_Mach_No"].asDouble();
    initialCond.u = InitialConditions["Initial_Conditions"]["V_1"].asDouble();
    initialCond.v = InitialConditions["Initial_Conditions"]["V_2"].asDouble();
    initialCond.T = InitialConditions["Initial_Conditions"]["Temperature_Static_Inlet"].asDouble();
    initialCond.T_inf = InitialConditions["Initial_Conditions"]["T_inf"].asDouble();
    initialCond.P_inf = InitialConditions["Initial_Conditions"]["P_inf"].asDouble();
    initialCond.Rho_inf = InitialConditions["Initial_Conditions"]["Rho_inf"].asDouble();

    return true;
}

// Reads initial and boundary conditions from a JSON file and populates the structures
bool readInitialAndBoundaryConditions(const std::string &filename, InitialCondition &initialCond, InletCondition &inletCond, ExitCondition &exitCond)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << "\n";
        return false;
    }

    Json::Reader reader;
    Json::Value data;
    if (!reader.parse(file, data))
    {
        throw std::runtime_error("Error: Failed to parse JSON input.");
    }

    const auto &IC = data["Initial_Conditions"];
    initialCond.P = IC["Pressure_Static_Inlet"].asDouble();
    initialCond.Rho = IC["Rho_Static_Inlet"].asDouble();
    initialCond.M = IC["Inlet_Mach_No"].asDouble();
    initialCond.u = IC["V_1"].asDouble();
    initialCond.v = IC["V_2"].asDouble();
    initialCond.T = IC["Temperature_Static_Inlet"].asDouble();
    initialCond.T_inf = IC["T_inf"].asDouble();
    initialCond.P_inf = IC["P_inf"].asDouble();
    initialCond.Rho_inf = IC["Rho_inf"].asDouble();

    const auto &Inlet = data["Boundary_Conditions"]["Inlet"];
    inletCond.P = Inlet["Pressure_Static_Inlet"].asDouble();
    inletCond.Rho = Inlet["Rho_Static_Inlet"].asDouble();
    inletCond.M = Inlet["Inlet_Mach_No"].asDouble();
    inletCond.u = Inlet["V_1"].asDouble();
    inletCond.v = Inlet["V_2"].asDouble();
    inletCond.T = Inlet["Temperature_Static_Inlet"].asDouble();
    inletCond.T_inf = Inlet["T_inf"].asDouble();
    inletCond.P_inf = Inlet["P_inf"].asDouble();
    inletCond.Rho_inf = Inlet["Rho_inf"].asDouble();

    const auto &Exit = data["Boundary_Conditions"]["Exit"];
    exitCond.P = Exit["exit_pressure"].asDouble();
    exitCond.Rho = Exit["exit_density"].asDouble();
    exitCond.T = Exit["exit_temperature"].asDouble();
    exitCond.M = Exit["exit_mach_no"].asDouble();
    exitCond.u = Exit["V_1"].asDouble();
    exitCond.v = Exit["V_2"].asDouble();

    return true;
}

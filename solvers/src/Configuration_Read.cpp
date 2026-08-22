#include "definitions.h"
#include "Globals.h"
#include "Directory_Files.h"
#include "Test_Cases.h"
#include "LES_Models.h"

#include <cctype>
#include <cmath>
#include <iostream>

/// @brief
/// @param jsonFileName  The JSON file name
/// @return void     No return value
/// @details This function reads the JSON file and sets the parameters for the solver and all other parameters
void readJSON(const std::string &jsonFileName)
{
    std::ifstream fileStream = CFD_OpenInputFileOrThrow(jsonFileName, "readJSON");

    Json::Reader reader;
    Json::Value root;
    if (!reader.parse(fileStream, root))
    {
        throw std::runtime_error("readJSON: failed to parse JSON input: " + jsonFileName);
    }
    if (!root.isMember("TestCase") || !root.isMember("Simulation") || !root.isMember("Solver"))
    {
        throw std::runtime_error("readJSON: missing required TestCase, Simulation, or Solver section in " + jsonFileName);
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
    if (Test_Case_JSON_File.empty())
        throw std::runtime_error("readJSON: TestCase.Test_Case_Json is empty");
    if (!std::filesystem::exists(Test_Case_JSON_File))
    {
        const std::filesystem::path configDir = std::filesystem::path(jsonFileName).parent_path();
        const std::filesystem::path candidate = configDir / Test_Case_JSON_File;
        if (!configDir.empty() && std::filesystem::exists(candidate))
            Test_Case_JSON_File = candidate.lexically_normal().string();
    }
    CFD_EnsureInputFileExists(Test_Case_JSON_File, "readJSON Test_Case_Json");

    // Simulation settings
    Initialize_Type = simData["Initialize_Type"].asInt();
    Is_Implicit_Method = simData["Is_Implicit_Method"].asBool();
    Total_Iterations = simData["Total_Iterations"].asInt();
    CFL = simData["CFL"].asDouble();
    CFD_RequireFinitePositive(CFL, "Simulation.CFL");
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
    Is_Char = solverData.isMember("Is_Char") ? solverData["Is_Char"].asBool() : false;
    Dissipation_Type = solverData["Dissipation_Type"].asInt();
    if (Total_Iterations < 0)
        throw std::runtime_error("Simulation.Total_Iterations must be non-negative");
    if (NUM_FLUX_COMPONENTS <= 0 || NUM_FLUX_COMPONENTS > 4)
        throw std::runtime_error("Solver.NUM_FLUX_COMPONENTS must be in [1, 4]");
    Is_MOVERS_1 = solverData["Is_MOVERS_1"].asBool();
    if (Flux_Type != Dissipation_Type)
    {
        std::cerr << "Note: Flux_Type (" << Flux_Type
                  << ") is only used for output directory names; Dissipation_Type (" << Dissipation_Type
                  << ") selects the numerical flux (see Net_Flux.cpp). Align them in JSON to avoid confusion.\n";
    }
    Enable_Entropy_Fix = solverData["Enable_Entropy_Fix"].asBool();
    if (solverData.isMember("Enable_AMR"))
        Enable_AMR = solverData["Enable_AMR"].asBool();
    if (solverData.isMember("AMR_Period"))
        AMR_Period = solverData["AMR_Period"].asInt();
    if (solverData.isMember("AMR_Start_Iteration"))
        AMR_Start_Iteration = solverData["AMR_Start_Iteration"].asInt();
    if (solverData.isMember("AMR_Gradient_Threshold"))
        AMR_Gradient_Threshold = solverData["AMR_Gradient_Threshold"].asDouble();
    if (solverData.isMember("AMR_Max_Fraction"))
        AMR_Max_Fraction = solverData["AMR_Max_Fraction"].asDouble();
    AMR_Coarsen_Threshold = 0.4 * AMR_Gradient_Threshold;
    if (solverData.isMember("AMR_Coarsen_Threshold"))
        AMR_Coarsen_Threshold = solverData["AMR_Coarsen_Threshold"].asDouble();
    // Current WENO reconstruction assumes structured uniform stencils.
    // AMR introduces hanging/coarse-fine interfaces and breaks that assumption.
    if (Is_WENO && Enable_AMR)
    {
        std::cerr << "Warning: WENO + AMR is not stencil-compatible in current implementation; "
                     "auto-disabling AMR to preserve shock structure.\n";
        Enable_AMR = false;
    }

    /* LES block is applied after Initialize_TestCase (needs cell arrays). */
    Store_Pending_LES_Config(root);

    // Limiter coefficients
    Limiter_Zeta = limiterData["Limiter_Zeta"].asDouble();
    Limiter_Zeta1 = limiterData["Limiter_Zeta1"].asDouble();
    /* Optional MUSCL damping (1 = no extra scaling). Lower toward 0 for stiffer problems at fixed CFL. */
    Second_Order_Limiter_Scale = 1.0;
    if (limiterData.isMember("Second_Order_Limiter_Scale"))
        Second_Order_Limiter_Scale = limiterData["Second_Order_Limiter_Scale"].asDouble();
    if (Second_Order_Limiter_Scale < 0.0)
        Second_Order_Limiter_Scale = 0.0;
    if (Second_Order_Limiter_Scale > 1.0)
        Second_Order_Limiter_Scale = 1.0;

    Second_Order_Phi_Blend = 1.0;
    if (limiterData.isMember("Second_Order_Phi_Blend"))
        Second_Order_Phi_Blend = limiterData["Second_Order_Phi_Blend"].asDouble();
    if (Second_Order_Phi_Blend < 0.0)
        Second_Order_Phi_Blend = 0.0;
    if (Second_Order_Phi_Blend > 1.0)
        Second_Order_Phi_Blend = 1.0;

    Second_Order_Dissipation_Blend = 1.0;
    if (limiterData.isMember("Second_Order_Dissipation_Blend"))
        Second_Order_Dissipation_Blend = limiterData["Second_Order_Dissipation_Blend"].asDouble();
    if (Second_Order_Dissipation_Blend < 0.0)
        Second_Order_Dissipation_Blend = 0.0;
    if (Second_Order_Dissipation_Blend > 1.0)
        Second_Order_Dissipation_Blend = 1.0;

    Venkat_K = 5.0;
    if (limiterData.isMember("Venkat_K"))
        Venkat_K = limiterData["Venkat_K"].asDouble();
    if (Venkat_K < 1e-6)
        Venkat_K = 1e-6;
    // Epsilon = limiterData["Epsilon"].asDouble(); // new
}

// Function to read parameters from a configuration file
map<string, string> ReadConfigFile(const string &filename)
{
    map<string, string> config;
    ifstream file = CFD_OpenInputFileOrThrow(filename, "ReadConfigFile");

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
    std::ifstream fileStream = CFD_OpenInputFileOrThrow(jsonFileName, "parseTestCaseJSON");

    Json::Reader reader;
    Json::Value root;
    if (!reader.parse(fileStream, root))
    {
        throw std::runtime_error("parseTestCaseJSON: failed to parse Test Case JSON: " + jsonFileName);
    }
    if (!root.isMember("GeneralDetails") || !root.isMember("mesh_parameters") || !root.isMember("Flow_Conditions"))
    {
        throw std::runtime_error("parseTestCaseJSON: missing required GeneralDetails, mesh_parameters, or Flow_Conditions section in " + jsonFileName);
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
    Grid_Size = meshParams.gridSize;
    if (meshParams.nx <= 1 || meshParams.ny <= 1 || Grid_Size <= 0)
        throw std::runtime_error("parseTestCaseJSON: mesh Nx, Ny, and Grid_Size must be positive and form at least one cell");

    // Flow Conditions
    const auto &flow = root["Flow_Conditions"];

    const auto &inlet = flow["inlet_conditions"];
    inletCond.type = inlet["inlet_type"].asString();
    inletCond.P = inlet["Pressure_Static_Inlet"].asDouble();
    inletCond.Rho = inlet["Rho_Static_Inlet"].asDouble();
    inletCond.M = inlet["Inlet_Mach_No"].asDouble();
    /* Accept either "u"/"v" or legacy "V_1"/"V_2" keys. */
    if (inlet.isMember("u") || inlet.isMember("v"))
    {
        inletCond.u = inlet.isMember("u") ? inlet["u"].asDouble() : 0.0;
        inletCond.v = inlet.isMember("v") ? inlet["v"].asDouble() : 0.0;
    }
    else
    {
        inletCond.u = inlet.isMember("V_1") ? inlet["V_1"].asDouble() : 0.0;
        inletCond.v = inlet.isMember("V_2") ? inlet["V_2"].asDouble() : 0.0;
    }
    /* Dimensional T uses R_GC; nondim T uses R=1/(gamma M^2) once M is known. */
    if (inletCond.Rho > 0.0 && inletCond.P > 0.0)
    {
        if (Non_Dimensional_Form && inletCond.M > 0.0)
            inletCond.T = inletCond.P / (inletCond.Rho * (1.0 / (gamma * inletCond.M * inletCond.M)));
        else
            inletCond.T = inletCond.P / (inletCond.Rho * R_GC);
    }
    /* If Mach is set but velocity is ~0, rebuild freestream speed from a=sqrt(gamma P/rho). */
    {
        const double a = (inletCond.Rho > 0.0)
                             ? std::sqrt(gamma * inletCond.P / inletCond.Rho)
                             : 0.0;
        const double speed = std::sqrt(inletCond.u * inletCond.u + inletCond.v * inletCond.v);
        if (inletCond.M > 0.5 && a > 0.0 && speed < 1.0e-12)
        {
            /* Half-cylinder freestream is -y; other cases often +x. Prefer -y when |v| key absent. */
            inletCond.u = 0.0;
            inletCond.v = -inletCond.M * a;
            std::cerr << "parseTestCaseJSON: rebuilt inlet velocity from Mach: u="
                      << inletCond.u << " v=" << inletCond.v << std::endl;
        }
    }

    const auto &exit = flow["exit_conditions"];
    exitCond.type = exit["exit_type"].asString();
    exitCond.P = exit["Pressure_Static_Exit"].asDouble();
    exitCond.Rho = exit["Rho_Static_Exit"].asDouble();
    exitCond.M = exit["Exit_Mach_No"].asDouble();
    exitCond.u = exit["u"].asDouble();
    exitCond.v = exit["v"].asDouble();

    const auto &wall = flow["wall_conditions"];
    wallCond.type = wall.isMember("wall_type") ? wall["wall_type"].asString() : "slip";
    wallCond.T = wall.isMember("wall_temperature") ? wall["wall_temperature"].asDouble() : 0.0;
    wallCond.u = wall.isMember("wall_velocity") ? wall["wall_velocity"].asDouble() : 0.0;
    wallCond.v = wall.isMember("wall_velocity") ? wall["wall_velocity"].asDouble() : 0.0;

    if (flow.isMember("Reynolds_number") && flow["Reynolds_number"].asDouble() > 0.0)
    {
        Re = flow["Reynolds_number"].asDouble();
        Inv_Re = 1.0 / Re;
    }
    if (flow.isMember("Prandtl_number") && flow["Prandtl_number"].asDouble() > 0.0)
    {
        Pr = flow["Prandtl_number"].asDouble();
        Inv_Pr = 1.0 / Pr;
    }

    const auto &init = flow["Initial_Conditions"];
    initCond.P = init["Pressure_Static_Inlet"].asDouble();
    initCond.Rho = init["Rho_Static_Inlet"].asDouble();
    initCond.M = init["Inlet_Mach_No"].asDouble();
    initCond.u = init["u"].asDouble();
    initCond.v = init["v"].asDouble();

    /* Route inlet/exit BC to supersonic vs subsonic implementations (Globals_Stub defaults are subsonic). */
    auto bcTypeIsSupersonic = [](const std::string &t) {
        std::string s = t;
        for (char &c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s == "supersonic";
    };
    Is_Inlet_SubSonic = !bcTypeIsSupersonic(inletCond.type);
    Is_Exit_SubSonic = !bcTypeIsSupersonic(exitCond.type);
    /* Driver JSON Is_Viscous selects NS vs Euler; keep wall BC flag aligned. */
    Is_Viscous_Wall = Is_Viscous;

    CFD_RequireFinitePositive(inletCond.P, "Flow_Conditions.inlet_conditions.Pressure_Static_Inlet");
    CFD_RequireFinitePositive(inletCond.Rho, "Flow_Conditions.inlet_conditions.Rho_Static_Inlet");
    CFD_RequireFinitePositive(exitCond.P, "Flow_Conditions.exit_conditions.Pressure_Static_Exit");
    CFD_RequireFinitePositive(exitCond.Rho, "Flow_Conditions.exit_conditions.Rho_Static_Exit");
    CFD_RequireFinitePositive(initCond.P, "Flow_Conditions.Initial_Conditions.Pressure_Static_Inlet");
    CFD_RequireFinitePositive(initCond.Rho, "Flow_Conditions.Initial_Conditions.Rho_Static_Inlet");
}

// Reads boundary conditions from a JSON file and populates the inlet and exit structures
bool readBoundaryConditionsjson(const std::string &filename, InletCondition &inletCond, ExitCondition &exitCond)
{
    std::ifstream file = CFD_OpenInputFileOrThrow(filename, "readBoundaryConditionsjson");
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
    std::ifstream file = CFD_OpenInputFileOrThrow(filename, "readInitialConditionsjson");
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
    std::ifstream file = CFD_OpenInputFileOrThrow(filename, "readInitialAndBoundaryConditions");

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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <iostream>

// ---------------------------------------------------------
// Defintions and Grid Readers
// ---------------------------------------------------------
#include "../Basic_Files/definitions.h"
#include "../Grid_Readers/Grid_Computations.cpp"
#include "../Grid_Readers/Co_Volume_Grid_Computations.cpp"

// ---------------------------------------------------------
// Solver Components
// ---------------------------------------------------------
#include "../Basic_Files/Initialize.cpp"
#include "../Basic_Files/Basic_Functions.cpp"
#include "../Basic_Files/Limiters.cpp"
#include "../Basic_Files/Error_Estimate_Update.cpp"
#include "../Basic_Files/Solver.cpp"
#include "../Basic_Files/output_files.cpp"
#include "../Basic_Files/Gradients.cpp"
#include "../Basic_Files/Utility_Files.cpp"
#include "../Basic_Files/WENO2D.cpp"
#include "../Basic_Files/Assemble_Matrix.cpp"

// ---------------------------------------------------------
// Flux Methods
// ---------------------------------------------------------
#include "../Flux_Types/Net_Flux.cpp"
#include "../Flux_Types/Average_Interface_Flux.cpp"
#include "../Flux_Types/Van_Leer.cpp"
#include "../Flux_Types/Dissipation.cpp"
#include "../Flux_Types/MOVERS_NWSC.cpp"
#include "../Flux_Types/LLF.cpp"
#include "../Flux_Types/RICCA.cpp"
#include "../Flux_Types/Roe_Scheme.cpp"
#include "../Flux_Types/Flux_Jacobian.cpp"

// ---------------------------------------------------------
// Grid Processing
// ---------------------------------------------------------
#include "../Grid_Readers/Create_Vtk_File.cpp"

// ---------------------------------------------------------
// Numerical Methods
// ---------------------------------------------------------
#include "../Numerical_Methods/Numerical_Method.cpp"

// ---------------------------------------------------------
// Boundary Conditions
// ---------------------------------------------------------
#include "../Boundary_Conditions/Inlet_Boundary_Conditions.cpp"
#include "../Boundary_Conditions/Exit_Boundary_Conditions.cpp"
#include "../Boundary_Conditions/Far_Field_Boundary_Condition.cpp"
#include "../Boundary_Conditions/Wall_Boundary_Conditions.cpp"
#include "../Boundary_Conditions/Boundary_Conditions.cpp"

// ---------------------------------------------------------
// Test Cases
// ---------------------------------------------------------
#include "../Test_Cases/Channel_Flow.cpp"
#include "../Test_Cases/SWBLI.cpp"
#include "../Test_Cases/Slip_Flow_2D.cpp"
#include "../Test_Cases/Shock_Bubble_Interaction.cpp"
#include "../Test_Cases/Shock_Tube_2D.cpp"
#include "../Test_Cases/Shock_Reflection_2D.cpp"
#include "../Test_Cases/Shock_Diffraction.cpp"
#include "../Test_Cases/Shock_Wedge_Reflection.cpp"
#include "../Test_Cases/Forward_Facing_Step.cpp"
#include "../Test_Cases/Ramp_15_Degree.cpp"
#include "../Test_Cases/Half_Cylinder_Test_Case.cpp"
#include "../Test_Cases/Flow_Over_Bump.cpp"
#include "../Test_Cases/Flat_Plate_Boundary_Layer.cpp"
#include "../Test_Cases/OED.cpp"
#include "../Test_Cases/Air_Foil.cpp"
#include "../Test_Cases/Scram_Jet_Inlet_Test_Case.cpp"
#include "../Test_Cases/Ellinga_TestCase.cpp"

// ---------------------------------------------------------
// Main Solver
// ---------------------------------------------------------
#include "../Euler_Solver/Main.cpp"
#include "../Euler_Solver/Configuration_Read.cpp"

#undef R // Avoid macro conflict with Pybind11

namespace py = pybind11;

using V_D = std::vector<double>;

PYBIND11_MODULE(CFD_Solver, m)
{
    m.doc() = "Python wrapper for CFD Solver";
   
    m.def("readJSON", &readJSON, "Reads the JSON input file and sets the global variables");

    m.def("testCase", 
      [](int &Test_Case) { testCase(Test_Case); }, 
      py::arg("Test_Case").noconvert(), 
      "Function to set test case ID");

    m.def("Set_DelU", &Set_DelU, py::arg("dU"), "Sets the delU values for the cells");

    // Geometric functions
    m.def("ReadGrid", &Read_Grid, "Reads the input grid file");
    m.def("FormCells", &Form_Cells, "Constructs cells from the grid data");
    m.def("ConstructCell", (void (*)(V_D &, V_D &, V_D &, V_D &))&Construct_Cell, "Constructs a cell");
    m.def("ConstructFace", &Construct_Face, "Evaluates face area, normals, and cell volume");
    m.def("ConstructCoVolumes", &Construct_Co_Volumes, "Constructs co-volumes for a cell");
    m.def("CheckCells", &Check_Cells, "Checks if the cells are valid and no negative volumes are present");
    m.def("ConstructGhostCells", &Construct_Ghost_Cells, "Construct ghost cells for the grid");
    m.def("WriteCellInfo", py::overload_cast<const std::string &>(&Write_Cell_Info), "Writes cell information to a file");

    m.def("Inviscid_Solver", &Inviscid_Solver,
          py::arg("error_file"), py::arg("solution_file"),
          "Solves using the inviscid solver with specified error and solution files.");

    m.def("Viscous_Solver", &Viscous_Solver,
          py::arg("error_file"), py::arg("solution_file"),
          "Solves using the viscous solver with specified error and solution files.");

    // For the methods that take no arguments
    m.def("Explicit_Method", &Explict_Method, "Executes the explicit method for solving");
    m.def("Runge_Kutta_Method", &Runge_Kutta_Method, "Executes the Runge-Kutta method for solving");

    // Time step functions (with arguments)
    m.def("Evaluate_Time_Step", [](int cellIndex)
          { Evaluate_Time_Step(cellIndex); }, "Evaluates the time step for the solver");
    m.def("Viscous_Time_Step_1", [](int cellIndex)
          { Viscous_Time_Step_1(cellIndex); }, "Viscous solver time step 1");
    m.def("Viscous_Time_Step_2", [](int cellIndex)
          { Viscous_Time_Step_2(cellIndex); }, "Viscous solver time step 2");
    m.def("Viscous_Time_Step_3", [](int cellIndex)
          { Viscous_Time_Step_3(cellIndex); }, "Viscous solver time step 3");
    m.def("Inviscid_Time_Step", [](int cellIndex)
          { Inviscid_Time_Step(cellIndex); }, "Inviscid solver time step");
    m.def("Directory_Name", &Directory_Name, "Adds Directory to file name");
    m.def("File_Name", &File_Name, "Adds Name to the file");

    m.def("Assemble_A", &Assemble_A, py::arg("A"), py::arg("dt"), py::return_value_policy::reference,
          "Assembles the global matrix A and vector b for CFD calculations");
    m.def("Assemble_A1", &Assemble_A1, py::arg("dt"), py::return_value_policy::reference,
          "Assembles only the non zero components of A and their associated indicies for CFD calculations");

    m.def("Assemble_b", &Assemble_b, py::arg("b"), py::return_value_policy::reference,
          "Assembles the global matrix A and vector b for CFD calculations");

    m.def("Compute_Flux_Jacobian", &Compute_Flux_Jacobian,
          py::arg("Cell_No"), py::arg("Ac"), py::arg("Face_No"),
          py::return_value_policy::reference,
          "Computes the flux Jacobian for a given face and cell");

    m.def("ComputeGhostCell_Flux_Jacobian", &ComputeGhostCell_Flux_Jacobian, py::arg("Ghost_Cell_No"),
          py::arg("Cell_No"), py::arg("Ac"), py::arg("Face_No"),
          py::return_value_policy::reference,
          "Computes the flux Jacobian for a given face and cell");

    m.def("get_NoPhysical_Cells", &get_NoPhysical_Cells, "Returns the number of physical cells");

    m.def("get_Min_dt", &get_Min_dt, "Returns the minimum time step of all the cells");
    m.def("get_Max_dt", &get_Max_dt, "Returns the maximum time step of all the cells");
    m.def("get_row_indices", &get_row_indices, "Returns list of row indicies");
    m.def("get_col_indices", &get_row_indices, "Returns list of col indicies");
    m.def("get_Values", &get_Values, "Returns list of values at that row and col");
    // Boundary conditions
    m.def("Apply_Boundary_Conditions", []()
          { Apply_Boundary_Conditions(); }, "Applies all boundary conditions based on the current setup");
    m.def("Estimate_Error", []()
          { Estimate_Error(); }, "Estimates the error in a given iteration");
    m.def("Update", []()
          { Update(); }, "Updates the Solution vector");

    m.def("ApplyBoundaryConditionsWithArgs", (void (*)(const int &, const bool &, const int &))&Apply_Boundary_Conditions, "Applies boundary conditions with arguments");
    m.def("InviscidWallBoundaryCondition", &InViscid_Wall_Boundary_Condition, "Applies inviscid wall boundary condition");
    m.def("ViscousWallBoundaryCondition", &Viscous_Wall_Boundary_Condition, "Applies viscous wall boundary condition");
    m.def("SymmetryBoundaryCondition", []()
          { Symmetry_Boundary_Condition(); }, "Applies symmetry boundary condition");
    m.def("SubsonicInletBoundaryCondition", &Subsonic_Inlet_Boundary_Condition, "Applies subsonic inlet boundary condition");
    m.def("SupersonicInletBoundaryCondition",
      py::overload_cast<>(&Supersonic_Inlet_Boundary_Condition),
      "Applies supersonic inlet boundary condition without parameters");

      m.def("SupersonicInletBoundaryConditionWithParams",
      py::overload_cast<double &, double &, double &, double &, double &>(&Supersonic_Inlet_Boundary_Condition),
      "Applies supersonic inlet boundary condition with specified parameters",
      py::arg("Pressure_Static_Inlet"), py::arg("Rho_Static_Inlet"),
      py::arg("u_Inlet"), py::arg("v_Inlet"), py::arg("Inlet_Mach_No"));
    // Flux and dissipation schemes
    m.def("EvaluateCellNetFlux1O", &Evaluate_Cell_Net_Flux_1O, "Evaluates the cell net flux using 1st order scheme");
    m.def("EvaluateCellNetFlux2O", &Evaluate_Cell_Net_Flux_2O, "Evaluates the cell net flux using 2nd order scheme");
    m.def("LaxFriedrichs", &Lax_Fedrichs, "Applies the Lax-Friedrichs method for dissipation");
    m.def("VanLeerFlux", &Van_Leer_Flux, "Applies Van Leer flux method");

    // Dissipation functions
    m.def("LLF", py::overload_cast<const int &, int &, const int &>(&LLF), "Applies LLF dissipation");
    m.def("ROE", py::overload_cast<const int &, int &, const int &>(&ROE), "Applies Roe dissipation");
    
    m.def("MOVERS", &MOVERS, "Applies MOVERS dissipation");

    // Gradient and reconstruction functions
    m.def("CalculatePrimitiveVariables", (void (*)(const int &, V_D &))&Calculate_Primitive_Variables, "Calculates primitive variables from U vector");
    m.def("CalculateComputationalVariables", (void (*)(V_D &))&Calculate_Computational_Variables, "Calculates computational variables from primitive variables");

    // Viscous flux functions
    m.def("Viscosity", &Viscosity, "Calculates the viscosity using Sutherland's law");
    m.def("ThermalConductivity", &Thermal_Conductivity, "Calculates thermal conductivity using Sutherland's law");
    m.def("EvaluateViscousFluxes", &Evaluate_Viscous_Fluxes, "Evaluates viscous fluxes using the Green-Gauss theorem");

    // Initialize and output functions
    m.def("Initialize", py::overload_cast<const int &>(&Initialize), "Initializes solver data");
    m.def("InitializeFromFile", py::overload_cast<const std::string &>(&Initialize), "Initializes solver from a file");

    m.def("WriteVTKFile", py::overload_cast<const std::string &>(&Write_VTK_File), "Writes VTK file for visualization");
    m.def("WriteCFFile", &Write_CF_File, "Writes the skin friction coefficient to a file");
    m.def("WriteErrorFile", &Write_Error_File, "Writes the error values to a file");

    // Bind overloaded Write_Solution functions
    m.def("Write_Solution_1", py::overload_cast<const std::string &, const int &>(&Write_Solution),
          py::arg("Op_file"), py::arg("type"), "Writes the solution data to a file based on the type.");

    // Read_Write_Grid and Append_Solution bindings
    m.def("Read_Write_Grid", &Read_Write_Grid, py::arg("Grid_Vtk_File"), py::arg("Final_Solution_File"),
          "Reads the grid file and appends the solution.");

    m.def("Append_Solution", &Append_Solution, py::arg("Solution_File"), py::arg("Final_Solution_File"),
          "Appends the solution to the final solution file.");

    // Bind overloaded Write_VTK_File functions
    m.def("Write_VTK_File_1", py::overload_cast<const std::string &, const std::string &>(&Write_VTK_File),
          py::arg("Op_file1"), py::arg("Op_file2"), "Writes VTK files for 3D data.");

    m.def("Write_VTK_File_2", py::overload_cast<const std::string &>(&Write_VTK_File),
          py::arg("File_Name"), "Writes VTK file for 2D data.");

    // Bind error, limiter, and CF file functions
    m.def("Write_Error_File", &Write_Error_File, py::arg("File_Name"), "Writes error data to a file.");
    m.def("Write_CF_File", &Write_CF_File, py::arg("File_Name"), "Writes CF data to a file.");

    m.def("Get_ErrorFileName", &Get_ErrorFileName, "Returns Error File name");
    m.def("Get_Initial_Solution_FileName", &Get_Initial_Solution_FileName, "Returns Initial Solution File");
    m.def("Get_Final_Solution_FileName", &Get_Final_Solution_FileName, "Returns Final Solution File");
    m.def("Get_Grid_VtkFile", &Get_Grid_VtkFile, "Returns GridFile File");
    m.def("Get_SolutionFile", &Get_SolutionFile, "Returns Initial Solution File");

    m.def("Write_Error_File", &Write_Error_File, "Writes Error to a file");
    m.def("Write_Solution", &Write_Solution, "Writes Solution to file");
    m.def("Read_Write_Grid", &Read_Write_Grid, "writes Grid_Vtk_File to Final_Solution_File");
    m.def("Append_Solution", &Append_Solution, "Writes Solution_File to Final_Solution_File");

    // Bind WriteMatrixToFile function
    m.def("WriteAMatrix", &Write_A_MatrixToFile, py::arg("A"), py::arg("File_Name"),
          "Writes a matrix and vector to a file.");
    // Bind WriteMatrixToFile function
    m.def("WriteMatrixToFile", &Write_b_VectorToFile, py::arg("b"), py::arg("File_Name"),
          "Writes a matrix and vector to a file.");

    // Wrapping of Test case file 
    m.def("Half_Cylinder_Flow", &Half_Cylinder_Flow, "Sets up and runs the half-cylinder flow simulation");
    m.def("Shock_Reflection_2d", &Shock_Reflection_2D, "Runs the shock reflection simulation in 2D");
    m.def("Ramp_15_Degree", &Ramp_15_Degree, "Simulates flow over a 15-degree ramp");
    m.def("Forward_Facing_Step", &Forward_Facing_Step, "Runs the forward-facing step simulation");
    m.def("Shock_Tube_2d", &Shock_Tube_2D, "Simulates the 2D shock tube case");
    m.def("Slip_Flow_2d", &Slip_Flow_2D, "Simulates slip flow in 2D");
    m.def("Channel_Flow", &Channel_Flow, "Simulates channel flow");
    m.def("Shock_Diffraction", &Shock_Diffraction, "Simulates shock diffraction");
    m.def("Shock_Wedge_Reflection", &Shock_Wedge_Reflection, "Simulates shock wedge reflection");
    m.def("Scram_Jet_Inlet", &Scram_Jet_Inlet, "Simulates scramjet inlet");
    m.def("Flow_Over_Bump", &Flow_Over_Bump, "Simulates flow over a bump");
    m.def("Flat_Plate_Boundary_Layer", &Flat_Plate_Boundary_Layer, "Simulates flow over a flat plate with boundary layer");
    m.def("SWBLI", &SWBLI, "Simulates shock wave boundary layer interaction (SWBLI)");
    m.def("Shock_Bubble_Interaction", &Shock_Bubble_Interaction, "Simulates shock bubble interaction");
    m.def("OED", &OED, "Simulates the OED test case");
    m.def("Ellinga_Carbuncle", &Ellinga_Carbuncle, "Simulates Ellinga carbuncle test case");

    // Exposing global variables in pybind11
    m.attr("Grid_Type") = &Grid_Type;
    m.attr("Initialize_Type") = &Initialize_Type;
    m.attr("Total_Iterations") = &Total_Iterations;
    m.attr("Limiter_Case") = &Limiter_Case;
    m.attr("Area_Weighted_Average") = &Area_Weighted_Average;
    m.attr("Flux_Type") = &Flux_Type;
    m.attr("Is_Second_Order") = &Is_Second_Order;
    m.attr("Time_Accurate") = &Time_Accurate;
    m.attr("Local_Time_Stepping") = &Local_Time_Stepping;
    m.attr("Non_Dimensional_Form") = &Non_Dimensional_Form;
    m.attr("Is_WENO") = &Is_WENO;
    m.attr("Is_Char") = &Is_Char;
    m.attr("Dissipation_Type") = &Dissipation_Type;
    m.attr("Is_MOVERS_1") = &Is_MOVERS_1;
    m.attr("Enable_Entropy_Fix") = &Enable_Entropy_Fix;
    m.attr("Test_Case") = &Test_Case;
    m.attr("Grid_Size") = &Grid_Size;
    m.attr("CFL") = &CFL;
    m.attr("Is_Viscous_Wall") = &Is_Viscous_Wall;
    //    m.attr("timer") = &timer;

    // Wrapping necessary global variables
    m.attr("Pressure_Static_Inlet") = &Pressure_Static_Inlet;
    m.attr("Inlet_Mach_No") = &Inlet_Mach_No;
    m.attr("Rho_Static_Inlet") = &Rho_Static_Inlet;
    m.attr("Temperature_Static_Inlet") = &Temperature_Static_Inlet;

    // Assuming the following variables are being used globally in the C++ file:
    // Wrap them so they can be accessed or modified from Python

    m.attr("Total_No_Cells") = py::cast(&Total_No_Cells);
    m.attr("No_Cartesian_Cells") = py::cast(&No_Cartesian_Cells);
    m.attr("No_Polar_Cells") = py::cast(&No_Polar_Cells);
    m.attr("No_Physical_Cells") = py::cast(&No_Physical_Cells);
    m.attr("No_Ghost_Cells") = py::cast(&No_Ghost_Cells);
    m.attr("Cells_in_Plane") = py::cast(&Cells_in_Plane);
    m.attr("Grid_Type") = py::cast(&Grid_Type);
    m.attr("global_temp") = py::cast(&global_temp);
    m.attr("R_Mid_dot_A") = py::cast(&R_Mid_dot_A);
    m.attr("Cell_Minimum_Length") = py::cast(&Cell_Minimum_Length);
    m.attr("CFL") = py::cast(&CFL);

    // Expose the string variables
    m.attr("Grid_File") = &Grid_File;
    m.attr("Initial_Solution_File") = &Initial_Solution_File;
    m.attr("Solution_File") = &Solution_File;
    m.attr("Error_File") = &Error_File;
    m.attr("Limiter_File") = &Limiter_File;
    m.attr("Final_Solution_File") = &Final_Solution_File;
    m.attr("Grid_Vtk_File") = &Grid_Vtk_File;
    m.attr("CF_File") = &CF_File;

    // Wrapping vectors and matrices defined in the C++ file for Python access
    m.attr("Cell_Face_Normals") = py::cast(&Cell_Face_Normals);
    m.attr("Cell_Face_Areas") = py::cast(&Cell_Face_Areas);
    m.attr("Cells_Diagonal_Vector") = py::cast(&Cells_Diagonal_Vector);
    m.attr("Cells_Cell_Center") = py::cast(&Cells_Cell_Center);
    m.attr("Cells_Center_Distances") = py::cast(&Cells_Center_Distances);
    m.attr("Distance_Bw_Cell_Centers") = py::cast(&Distance_Bw_Cell_Centers);
    m.attr("Cells_Vertices") = py::cast(&Cells_Vertices);
    m.attr("Cells_Area") = py::cast(&Cells_Area);
    m.attr("Cells_Inv_Area") = py::cast(&Cells_Inv_Area);
    m.attr("Face_Area_Components") = py::cast(&Face_Area_Components);
    m.attr("Face_Normal_Components") = py::cast(&Face_Normal_Components);
    m.attr("Global_Vec") = py::cast(&Global_Vec);
    m.attr("Cells_Area_Mag") = py::cast(&Cells_Area_Mag);
    m.attr("Cell_Face_Distances") = py::cast(&Cell_Face_Distances);
    m.attr("Vertices") = py::cast(&Vertices);
    m.attr("Cells_Volume") = py::cast(&Cells_Volume);
    m.attr("Cells_Inv_Volume") = py::cast(&Cells_Inv_Volume);
    m.attr("Distance_Ratio") = py::cast(&Distance_Ratio);

    // Wrapping integer lists like Cell Neighbours
    m.attr("Cell_Neighbour_indices") = py::cast(&Cell_Neighbour_indicies);
    m.attr("Cell_Neighbours") = py::cast(&Cell_Neighbours);
    m.attr("Cells_Secondary_Neighbours") = py::cast(&Cells_Secondary_Neighbours);

    // Wrapping boundary conditions lists
    m.attr("Wall_Cells_List") = py::cast(&Wall_Cells_List);
    m.attr("Inlet_Cells_List") = py::cast(&Inlet_Cells_List);
    m.attr("Exit_Cells_List") = py::cast(&Exit_Cells_List);
    m.attr("Symmetry_Cells_List") = py::cast(&Symmetry_Cells_List);
    m.attr("List_Dissipative_Cells") = py::cast(&List_Dissipative_Cells);
    m.attr("List_Non_Dissipative_Cells") = py::cast(&List_Non_Dissipative_Cells);
    m.attr("Far_Field_Out_Flow_List") = py::cast(&Far_Field_Out_Flow_List);
    m.attr("Far_Field_InFlow_List") = py::cast(&Far_Field_InFlow_List);

    // Wrapping boolean flags for the grid configurations
    m.attr("Is_Viscous_Wall") = py::cast(&Is_Viscous_Wall);
    m.attr("Is_2D_Flow") = py::cast(&Is_2D_Flow);
    m.attr("Is_Inlet_SubSonic") = py::cast(&Is_Inlet_SubSonic);
    m.attr("Is_Exit_SubSonic") = py::cast(&Is_Exit_SubSonic);
    m.attr("Enable_Far_Field") = py::cast(&Enable_Far_Field);
    m.attr("has_Symmetry_BC") = py::cast(&has_Symmetry_BC);
}

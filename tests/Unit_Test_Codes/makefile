# Executable name and directory
EXEC := Test_Gmsh_Reader
EXEC_DIRECTORY := ./
EXECUTABLES := Test_Gmsh_Reader

# Paths
BASE_PATH := ../
Solver_PATH := Basic_Files/
BC_PATH := Boundary_Conditions/
Flux_PATH:= Flux_Types/
Grid_PATH:= Grid_Readers/
Numerical_Methods_PATH:= Numerical_Methods/
Output_Files_PATH:= Output_Files/
Test_Cases_PATH:= Test_Cases/
Utility_Files_PATH:= Utility_Files/
UnitTests_PATH:= Unit_Test_Codes/
MAIN_PATH := Euler_Solver/
IPATH := MYLIBRARY/include/
SPATH := MYLIBRARY/src/
LPATH := ../

# Compilation flags
IFLAGS := -I $(IPATH) -I $(LPATH)
IFLAGS += -I /opt/homebrew/include
IFLAGS += -I /opt/homebrew/include/vtk-9.4
IFLAGS += -I /opt/homebrew/Cellar/jsoncpp/1.9.6/include
IFLAGS += -I/opt/homebrew/include/json

#LFLAGS := -L$(LPATH)  # Removed unnecessary -lnlohmann_json
LFLAGS := -L/opt/homebrew/lib -ljsoncpp
LFLAGS += -L/opt/homebrew/lib \
          -lvtkCommonCore-9.4 \
          -lvtkCommonDataModel-9.4 \
          -lvtkCommonExecutionModel-9.4 \
          -lvtkFiltersCore-9.4 \
          -lvtkFiltersGeometry-9.4 \
          -lvtkIOXML-9.4 \
          -lvtkRenderingCore-9.4 \
          -lvtkIOCore-9.4 \
          -lvtkIOImage-9.4 \
          -lvtkInteractionStyle-9.4 \
          -lvtkRenderingOpenGL2-9.4 \
          -lvtksys-9.4
DFLAGS := -ggdb3 -g -std=c++17
OFLAGS := -O3 -O2 -O
WFLAGS := -Wall -Wno-long-long -Wpedantic -Wextra
CFLAGS := -c $(IFLAGS) $(DFLAGS) $(OFLAGS) $(WFLAGS)

# Compiler and tools
CC := g++
STRIP := strip

# Object files
# Object files
OBJS = $(BASE_PATH)$(Grid_PATH)Grid_Computations.o \
       $(BASE_PATH)$(Grid_PATH)Co_Volume_Grid_Computations.o \
       $(BASE_PATH)$(Grid_PATH)Read_Gmsh_File.o \
       $(BASE_PATH)$(Grid_PATH)Grid_Refine_Functions.o \
       $(BASE_PATH)$(Solver_PATH)Initialize.o \
       $(BASE_PATH)$(Solver_PATH)Basic_Functions.o \
       $(BASE_PATH)$(Solver_PATH)Limiters.o \
       $(BASE_PATH)$(Solver_PATH)Error_Estimate_Update.o \
       $(BASE_PATH)$(Solver_PATH)Solver.o \
       $(BASE_PATH)$(Solver_PATH)output_files.o \
       $(BASE_PATH)$(Solver_PATH)Utility_Files.o \
       $(BASE_PATH)$(Solver_PATH)WENO2D.o \
       $(BASE_PATH)$(Solver_PATH)Dir_Files.o \
       $(BASE_PATH)$(Flux_PATH)Net_Flux.o \
       $(BASE_PATH)$(Flux_PATH)Average_Interface_Flux.o \
       $(BASE_PATH)$(Flux_PATH)Van_Leer.o \
       $(BASE_PATH)$(Flux_PATH)MOVERS.o \
       $(BASE_PATH)$(Flux_PATH)MOVERS_NWSC.o \
       $(BASE_PATH)$(Flux_PATH)LLF.o \
       $(BASE_PATH)$(Flux_PATH)RICCA.o \
       $(BASE_PATH)$(Flux_PATH)Roe_Scheme.o \
       $(BASE_PATH)$(Flux_PATH)Flux_Jacobian.o \
       $(BASE_PATH)$(Grid_PATH)Create_Vtk_File.o \
       $(BASE_PATH)$(Numerical_Methods_PATH)Numerical_Method.o \
       $(BASE_PATH)$(BC_PATH)Inlet_Boundary_Conditions.o \
       $(BASE_PATH)$(BC_PATH)Exit_Boundary_Conditions.o \
       $(BASE_PATH)$(BC_PATH)Far_Field_Boundary_Condition.o \
       $(BASE_PATH)$(BC_PATH)Wall_Boundary_Conditions.o \
       $(BASE_PATH)$(BC_PATH)Boundary_Conditions.o \
       $(BASE_PATH)$(Solver_PATH)Assemble_Matrix.o \
       $(BASE_PATH)$(Test_Cases_PATH)Channel_Flow.o \
       $(BASE_PATH)$(Test_Cases_PATH)SWBLI.o \
       $(BASE_PATH)$(Test_Cases_PATH)Slip_Flow_2D.o \
       $(BASE_PATH)$(Test_Cases_PATH)Shock_Bubble_Interaction.o \
       $(BASE_PATH)$(Test_Cases_PATH)Shock_Tube_2D.o \
       $(BASE_PATH)$(Test_Cases_PATH)Shock_Reflection_2D.o \
       $(BASE_PATH)$(Test_Cases_PATH)Shock_Diffraction.o \
       $(BASE_PATH)$(Test_Cases_PATH)Shock_Wedge_Reflection.o \
       $(BASE_PATH)$(Test_Cases_PATH)Forward_Facing_Step.o \
       $(BASE_PATH)$(Test_Cases_PATH)Ramp_15_Degree.o \
       $(BASE_PATH)$(Test_Cases_PATH)Half_Cylinder_Test_Case.o \
       $(BASE_PATH)$(Test_Cases_PATH)Flow_Over_Bump.o \
       $(BASE_PATH)$(Test_Cases_PATH)Flat_Plate_Boundary_Layer.o \
       $(BASE_PATH)$(Test_Cases_PATH)OED.o \
       $(BASE_PATH)$(Test_Cases_PATH)Air_Foil.o \
       $(BASE_PATH)$(Test_Cases_PATH)Scram_Jet_Inlet_Test_Case.o \
       $(BASE_PATH)$(Test_Cases_PATH)Ellinga_TestCase.o \
       $(BASE_PATH)$(MAIN_PATH)Configuration_Read.o \
       $(BASE_PATH)$(UnitTests_PATH)Test_Gmsh_Reader.o

# Default target to build the executable
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o $(EXEC_DIRECTORY)$@

# Pattern rule for .o files from .cpp
%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJS)
	rm -f $(EXEC_DIRECTORY)$(EXEC)
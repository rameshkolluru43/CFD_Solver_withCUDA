# VTK Grid Tests

This folder contains all the VTK grid reading tests and related files developed for validating the CFD solver's grid processing capabilities.

## 📁 Contents

### Main Integration Test
- **`test_read_vtk_grid_integrated.cpp`** - Complete integration test with real CFD solver functions
- **`test_read_vtk_grid_integrated`** - Compiled executable

### Development Test Files
- **`test_read_vtk_grid.cpp`** - Original test file
- **`test_read_vtk_grid`** - Original compiled executable
- **`comprehensive_vtk_test.cpp`** - Comprehensive VTK testing
- **`full_vtk_test.cpp`** - Full VTK functionality test
- **`full_vtk_test`** - Compiled executable
- **`simple_vtk_test.cpp`** - Simple VTK reading test
- **`simple_vtk_test`** - Compiled executable
- **`simple_grid_test.cpp`** - Simple grid processing test
- **`simple_grid_test`** - Compiled executable
- **`test_vtk_reader.cpp`** - VTK reader testing
- **`test_load_mesh.cpp`** - Mesh loading test
- **`boolean_test.cpp`** - Boolean logic testing
- **`boolean_test`** - Compiled executable

### Build Files
- **`CMakeLists_test.txt`** - CMake configuration for tests
- **`Makefile.test`** - Makefile for test compilation

### Configuration Files
- **`VTK_Test_Config.json`** - Test configuration
- **`VTK_Test_Ramp.json`** - Ramp test case configuration

## 🚀 Running Tests

### Main Integration Test
```bash
# Compile
g++ -std=c++17 -o test_read_vtk_grid_integrated test_read_vtk_grid_integrated.cpp

# Run
./test_read_vtk_grid_integrated
```

### Other Tests
```bash
# Compile any test
g++ -std=c++17 -o <test_name> <test_name>.cpp

# Run
./<test_name>
```

## 📊 Test Results

The main integration test (`test_read_vtk_grid_integrated.cpp`) validates:

### ✅ Core Functionality
- VTK file reading (867 points, 800 cells)
- Real geometric calculations
- Face normal computation
- Area calculations
- Cell center computation

### ✅ Advanced Features
- Boundary data structures
- Ghost cell construction (131 cells)
- Co-Volume implementation
- Cell-Ghost connectivity
- Data structure validation

### 📈 Quality Metrics
- Grid quality: Area ratio 1.86 (well-conditioned)
- Boundary distribution: 16 inlet, 16 exit, 99 wall
- All geometric validations: PASSED

## 📚 Documentation

For complete test documentation, see:
- **`../TEST_SUMMARY.md`** - Comprehensive test results and analysis
- **`../README.md`** - Main project documentation with validation section

## 🔧 Development Notes

These tests evolved through 5 phases:
1. Basic VTK file reading
2. Mock implementation testing
3. Real CFD solver integration
4. Enhanced validation with Check_Cells()
5. Boundary data structures and ghost cells

The final integration test provides a complete validation framework for the CFD solver's grid processing pipeline.

---
*Test files organized: September 27, 2025*
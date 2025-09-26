#!/bin/bash

# Set Doxygen config file
DOXYFILE="Doxyfile_Cleaned"

# Define directories containing C++ and Python files
SRC_DIRS=(
    "Basic_Files"
    "Boundary_Conditions"
    "Euler_Solver"
    "Flux_Types"
    "Grid_Readers"
    "Numerical_Methods"
    "Test_Cases"
    "Python_Wrapper"
)

# Ensure Doxygen processes the correct files
echo "Updating Doxyfile with correct INPUT paths..."
sed -i '' "/^INPUT /c\INPUT = ${SRC_DIRS[*]}" "$DOXYFILE"
sed -i '' "/^FILE_PATTERNS /c\FILE_PATTERNS = *.cpp *.h *.hpp *.py" "$DOXYFILE"
sed -i '' "/^RECURSIVE /c\RECURSIVE = YES" "$DOXYFILE"

# Run Doxygen
echo "Running Doxygen..."
doxygen "$DOXYFILE"

# Check if documentation was generated
if [ -d "html" ]; then
    echo "Documentation successfully generated in html/index.html"
else
    echo "Error: Documentation was not generated. Check Doxygen settings."
fi
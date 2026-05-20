#!/bin/bash

# Set Doxygen config file
DOXYFILE="docs/Doxyfile_Cleaned"
LOG_FILE="docs/logs/doxygen_log.txt"

# Define directories containing C++ and CUDA files
SRC_DIRS=(
    "src"
    "include"
    "CUDA_KERNELS"
    "Basic_Function_Files"
    "docs"
    "Python_Wrapper"
)

echo "=== CFD Solver with Enhanced Flux Schemes Documentation Update ==="
echo "Updating documentation with production-ready flux schemes and CUDA kernels..."
echo "Enhanced Features: Van Leer, Roe (1st/2nd order), AUSM flux schemes"
echo "Source directories: ${SRC_DIRS[*]}"

# Ensure Doxygen processes the correct files with CUDA support
echo "Using pre-configured Doxyfile with CUDA support..."
# The Doxyfile_Cleaned is already configured with:
# - INPUT paths including CUDA_KERNELS directory
# - FILE_PATTERNS including *.cu and *.cuh files
# - EXTENSION_MAPPING for CUDA files
# - RECURSIVE enabled for directory traversal

# Run Doxygen
echo "Running Doxygen to generate documentation..."
doxygen "$DOXYFILE" 2>&1 | tee "$LOG_FILE"

# Check if documentation was generated
if [ -d "docs/doxygen/html" ]; then
    echo "✅ Documentation successfully generated!"
    echo "📁 Main documentation: docs/doxygen/html/index.html"
    echo "📊 Doxygen log saved to: $LOG_FILE"
    
    # Count files processed
    if [ -f "$LOG_FILE" ]; then
        file_count=$(grep -c "Parsing file" "$LOG_FILE" || echo "0")
        echo "📄 Files processed: $file_count"
        
        # Check if CUDA files were processed
        cuda_count=$(grep -c "\.cu\|\.cuh" "$LOG_FILE" 2>/dev/null || echo "0")
        if [ "${cuda_count}" -gt 0 ] 2>/dev/null; then
            echo "🚀 CUDA files processed: $cuda_count"
        fi
    fi
    
    echo ""
    echo "=== Enhanced Documentation Features ==="
    echo "• 🚀 Production-ready flux schemes (Van Leer, Roe, AUSM)"
    echo "• 🔧 Complete CFD solver API documentation"
    echo "• 💻 CUDA GPU acceleration kernels"
    echo "• 📊 Performance benchmarks and optimization guides"
    echo "• 🧮 Mathematical formulations and algorithms"
    echo "• ✅ Error handling and validation frameworks"
    echo "• 📚 Comprehensive technical implementation guides"
    echo "• 🎯 Integration examples and usage tutorials"
    echo ""
    echo "🌐 Open docs/doxygen/html/index.html in your web browser to view the documentation"
    
else
    echo "❌ Error: Documentation was not generated"
    echo "🔍 Check $LOG_FILE for detailed error messages"
    if [ -f "$LOG_FILE" ]; then
        echo "Last few lines of log:"
        tail -10 "$LOG_FILE"
    fi
fi
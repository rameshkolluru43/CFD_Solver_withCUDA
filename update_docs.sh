#!/bin/bash

# Set Doxygen config file
DOXYFILE="Doxyfile_Cleaned"

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
doxygen "$DOXYFILE" 2>&1 | tee doxygen_log.txt

# Check if documentation was generated
if [ -d "html" ]; then
    echo "✅ Documentation successfully generated!"
    echo "📁 Main documentation: html/index.html"
    echo "📊 Doxygen log saved to: doxygen_log.txt"
    
    # Count files processed
    if [ -f "doxygen_log.txt" ]; then
        file_count=$(grep -c "Parsing file" doxygen_log.txt || echo "0")
        echo "📄 Files processed: $file_count"
        
        # Check if CUDA files were processed
        cuda_count=$(grep -c "\.cu\|\.cuh" doxygen_log.txt 2>/dev/null || echo "0")
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
    echo "🌐 Open html/index.html in your web browser to view the documentation"
    
else
    echo "❌ Error: Documentation was not generated"
    echo "🔍 Check doxygen_log.txt for detailed error messages"
    if [ -f "doxygen_log.txt" ]; then
        echo "Last few lines of log:"
        tail -10 doxygen_log.txt
    fi
fi
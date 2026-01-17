#!/bin/bash

# CFD Solver Documentation Viewer
# This script opens the generated Doxygen documentation in the default browser

echo "=== CFD Solver Suite Documentation ==="
echo "Opening documentation in your default browser..."

# Get the absolute path to the HTML documentation
DOC_PATH="$(pwd)/html/index.html"

if [ -f "$DOC_PATH" ]; then
    echo "Documentation path: $DOC_PATH"
    
    # Open in default browser (macOS)
    if [[ "$OSTYPE" == "darwin"* ]]; then
        open "$DOC_PATH"
    # Open in default browser (Linux)
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        xdg-open "$DOC_PATH"
    # Open in default browser (Windows/WSL)
    elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
        start "$DOC_PATH"
    else
        echo "Please open the following file in your browser:"
        echo "file://$DOC_PATH"
    fi
    
    echo ""
    echo "=== Documentation Features ==="
    echo "✓ Search functionality enabled"
    echo "✓ Class hierarchy with inheritance graphs"
    echo "✓ File dependency diagrams"
    echo "✓ Function and variable cross-references"
    echo "✓ Source code browsing with syntax highlighting"
    echo "✓ Detailed API documentation for all components"
    echo ""
    echo "=== Key Sections to Explore ==="
    echo "• Main Page: Complete project overview and features"
    echo "• Classes: All C++ classes with detailed member documentation"
    echo "• Files: Source file organization and dependencies"
    echo "• Functions: Complete function reference with parameters"
    echo "• Search: Find any symbol, function, or documentation"
    
else
    echo "Error: Documentation not found at $DOC_PATH"
    echo "Please run 'doxygen Doxyfile_Cleaned' to generate documentation first."
    exit 1
fi
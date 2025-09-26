#!/bin/bash
# Build documentation script for CFD Solver with CUDA
# Author: Ramesh Kolluru
# Version: 0.1.0-alpha

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

echo -e "${BLUE}CFD Solver with CUDA - Documentation Builder${NC}"
echo "=============================================="

# Change to project root
cd "$PROJECT_ROOT"

# Check if Doxygen is installed
if ! command -v doxygen &> /dev/null; then
    echo -e "${RED}Error: Doxygen is not installed${NC}"
    echo "Please install Doxygen:"
    echo "  Ubuntu/Debian: sudo apt install doxygen"
    echo "  CentOS/RHEL:   sudo yum install doxygen"
    echo "  macOS:         brew install doxygen"
    exit 1
fi

echo -e "${GREEN}✓ Doxygen found: $(doxygen --version)${NC}"

# Check if dot (Graphviz) is available for diagrams
if command -v dot &> /dev/null; then
    echo -e "${GREEN}✓ Graphviz found: $(dot -V 2>&1 | head -n1)${NC}"
    HAVE_DOT="YES"
else
    echo -e "${YELLOW}⚠ Graphviz not found - diagrams will not be generated${NC}"
    echo "  To enable diagrams, install Graphviz:"
    echo "  Ubuntu/Debian: sudo apt install graphviz"
    echo "  CentOS/RHEL:   sudo yum install graphviz"
    echo "  macOS:         brew install graphviz"
    HAVE_DOT="NO"
fi

# Create temporary Doxyfile with updated HAVE_DOT setting
TEMP_DOXYFILE="$(mktemp)"
cp Doxyfile "$TEMP_DOXYFILE"
sed -i.bak "s/HAVE_DOT.*=.*/HAVE_DOT = $HAVE_DOT/" "$TEMP_DOXYFILE"

# Clean previous documentation
echo -e "${YELLOW}Cleaning previous documentation...${NC}"
rm -rf docs/html docs/latex

# Create docs directory if it doesn't exist
mkdir -p docs

# Generate documentation
echo -e "${YELLOW}Generating documentation...${NC}"
if doxygen "$TEMP_DOXYFILE"; then
    echo -e "${GREEN}✓ Documentation generated successfully${NC}"
else
    echo -e "${RED}✗ Documentation generation failed${NC}"
    rm -f "$TEMP_DOXYFILE" "$TEMP_DOXYFILE.bak"
    exit 1
fi

# Cleanup
rm -f "$TEMP_DOXYFILE" "$TEMP_DOXYFILE.bak"

# Check if documentation was created
if [ -d "docs/html" ] && [ -f "docs/html/index.html" ]; then
    echo -e "${GREEN}✓ HTML documentation created in docs/html/${NC}"
    
    # Count files generated
    HTML_FILES=$(find docs/html -name "*.html" | wc -l)
    echo "  Generated $HTML_FILES HTML files"
    
    # Check main page
    if [ -f "docs/html/index.html" ]; then
        echo "  Main page: docs/html/index.html"
    fi
    
    # Check for search functionality
    if [ -f "docs/html/search/search.js" ]; then
        echo -e "${GREEN}  ✓ Search functionality enabled${NC}"
    fi
    
    # Check for class diagrams
    if [ -d "docs/html" ] && find docs/html -name "*.png" | head -1 > /dev/null; then
        echo -e "${GREEN}  ✓ Class diagrams generated${NC}"
    fi
    
else
    echo -e "${RED}✗ HTML documentation not found${NC}"
    exit 1
fi

# Optional: Generate LaTeX documentation
if command -v pdflatex &> /dev/null && [ "$1" = "--pdf" ]; then
    echo -e "${YELLOW}Generating PDF documentation...${NC}"
    
    # Update Doxyfile for LaTeX generation
    TEMP_DOXYFILE="$(mktemp)"
    cp Doxyfile "$TEMP_DOXYFILE"
    sed -i.bak "s/GENERATE_LATEX.*=.*/GENERATE_LATEX = YES/" "$TEMP_DOXYFILE"
    
    if doxygen "$TEMP_DOXYFILE"; then
        cd docs/latex
        if make; then
            echo -e "${GREEN}✓ PDF documentation generated: docs/latex/refman.pdf${NC}"
        else
            echo -e "${YELLOW}⚠ PDF generation failed${NC}"
        fi
        cd "$PROJECT_ROOT"
    fi
    
    rm -f "$TEMP_DOXYFILE" "$TEMP_DOXYFILE.bak"
fi

# Generate statistics
echo ""
echo "Documentation Statistics:"
echo "========================"

if [ -d "docs/html" ]; then
    HTML_SIZE=$(du -sh docs/html | cut -f1)
    echo "HTML documentation size: $HTML_SIZE"
fi

if [ -f "docs/latex/refman.pdf" ]; then
    PDF_SIZE=$(du -sh docs/latex/refman.pdf | cut -f1)
    echo "PDF documentation size: $PDF_SIZE"
fi

# Count documented items
if [ -f "docs/html/functions.html" ]; then
    FUNCTIONS=$(grep -c "class=\"memItemLeft\"" docs/html/functions.html 2>/dev/null || echo "Unknown")
    echo "Documented functions: $FUNCTIONS"
fi

if [ -f "docs/html/annotated.html" ]; then
    CLASSES=$(grep -c "class=\"memItemLeft\"" docs/html/annotated.html 2>/dev/null || echo "Unknown")
    echo "Documented classes: $CLASSES"
fi

echo ""
echo -e "${GREEN}Documentation build complete!${NC}"
echo ""
echo "To view the documentation:"
echo "  Open: docs/html/index.html"
echo "  Or run: python3 -m http.server 8000 -d docs/html"
echo "  Then visit: http://localhost:8000"

# Optional: Open documentation in browser
if [ "$1" = "--open" ] || [ "$2" = "--open" ]; then
    if command -v xdg-open &> /dev/null; then
        xdg-open "docs/html/index.html"
    elif command -v open &> /dev/null; then
        open "docs/html/index.html"
    else
        echo "Please open docs/html/index.html manually"
    fi
fi
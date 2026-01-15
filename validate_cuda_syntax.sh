#!/bin/bash
# Syntax Validation Script for New CUDA Kernels
# Date: 2026-01-15
# Usage: ./validate_cuda_syntax.sh

echo "=========================================="
echo "CUDA Kernel Syntax Validation"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if CUDA is available
if ! command -v nvcc &> /dev/null; then
    echo -e "${RED}ERROR: nvcc not found${NC}"
    echo "CUDA Toolkit must be installed to validate syntax"
    echo "Please install CUDA Toolkit from: https://developer.nvidia.com/cuda-downloads"
    exit 1
fi

echo -e "${GREEN}✓ CUDA Toolkit found${NC}"
nvcc --version | head -1
echo ""

# Set paths
CUDA_DIR="CUDA_KERNELS"
BUILD_DIR="build/syntax_check"

# Create temporary build directory
mkdir -p $BUILD_DIR
cd $BUILD_DIR

echo "=========================================="
echo "Validating New Kernel Files"
echo "=========================================="
echo ""

# List of new files to validate
NEW_FILES=(
    "Roe_Flux_Cuda_Kernels.cu"
    "HLLC_LLF_Flux_Cuda_Kernels.cu"
    "MUSCL_WENO_Reconstruction_Cuda_Kernels.cu"
    "Advanced_Flux_Schemes_Cuda_Wrappers.cu"
    "Reconstruction_Schemes_Cuda_Wrappers.cu"
)

HEADERS=(
    "Advanced_Flux_Schemes_Cuda.h"
    "Reconstruction_Schemes_Cuda.h"
)

# Validation flags
PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0

# Validation function
validate_file() {
    local file=$1
    local filepath="../../$CUDA_DIR/$file"
    
    echo -n "Validating $file ... "
    
    # Check if file exists
    if [ ! -f "$filepath" ]; then
        echo -e "${RED}FAIL${NC} (file not found)"
        ((FAIL_COUNT++))
        return 1
    fi
    
    # Syntax check with nvcc (compile to PTX only, no linking)
    nvcc -std=c++17 -arch=sm_60 -ptx \
         -I../../$CUDA_DIR \
         -I../../include \
         --expt-relaxed-constexpr \
         -Werror all-warnings \
         "$filepath" -o "${file}.ptx" 2>&1 | tee "${file}.log"
    
    local exit_code=${PIPESTATUS[0]}
    
    if [ $exit_code -eq 0 ]; then
        # Check for warnings
        if grep -q "warning" "${file}.log"; then
            echo -e "${YELLOW}WARN${NC} (compiled with warnings)"
            echo "  See: $BUILD_DIR/${file}.log"
            ((WARN_COUNT++))
        else
            echo -e "${GREEN}PASS${NC}"
            ((PASS_COUNT++))
        fi
        
        # Report register usage
        if [ -f "${file}.ptx" ]; then
            local max_regs=$(grep -o "\.maxnreg [0-9]*" "${file}.ptx" | awk '{print $2}' | sort -n | tail -1)
            if [ ! -z "$max_regs" ]; then
                echo "  Max registers: $max_regs"
                if [ $max_regs -gt 64 ]; then
                    echo -e "  ${YELLOW}Warning: High register usage may reduce occupancy${NC}"
                fi
            fi
        fi
        
        return 0
    else
        echo -e "${RED}FAIL${NC}"
        echo "  Errors found. See: $BUILD_DIR/${file}.log"
        ((FAIL_COUNT++))
        
        # Show first error
        echo "  First error:"
        grep "error:" "${file}.log" | head -1 | sed 's/^/    /'
        
        return 1
    fi
}

# Validate headers
echo "Checking Headers..."
for header in "${HEADERS[@]}"; do
    filepath="../../$CUDA_DIR/$header"
    if [ -f "$filepath" ]; then
        echo -e "  $header: ${GREEN}EXISTS${NC}"
        
        # Check for common issues
        if ! grep -q "#ifndef" "$filepath"; then
            echo -e "    ${YELLOW}Warning: Missing include guard${NC}"
        fi
        if ! grep -q "cuda_runtime.h" "$filepath"; then
            echo -e "    ${YELLOW}Warning: Missing cuda_runtime.h${NC}"
        fi
    else
        echo -e "  $header: ${RED}MISSING${NC}"
        ((FAIL_COUNT++))
    fi
done
echo ""

# Validate implementation files
echo "Checking Implementation Files..."
for file in "${NEW_FILES[@]}"; do
    validate_file "$file"
    echo ""
done

echo "=========================================="
echo "Validation Summary"
echo "=========================================="
echo -e "Passed:   ${GREEN}$PASS_COUNT${NC}"
echo -e "Warnings: ${YELLOW}$WARN_COUNT${NC}"
echo -e "Failed:   ${RED}$FAIL_COUNT${NC}"
echo ""

# Check CMakeLists.txt
echo "Checking CMakeLists.txt..."
if grep -q "Roe_Flux_Cuda_Kernels.cu" "../../CMakeLists.txt"; then
    echo -e "  ${GREEN}✓ New files added to CMakeLists.txt${NC}"
else
    echo -e "  ${RED}✗ New files NOT in CMakeLists.txt${NC}"
    echo "  Run: vi ../../CMakeLists.txt"
    ((FAIL_COUNT++))
fi
echo ""

# Overall result
echo "=========================================="
if [ $FAIL_COUNT -eq 0 ]; then
    echo -e "${GREEN}✓ ALL VALIDATIONS PASSED${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. cd ../../build"
    echo "  2. cmake .."
    echo "  3. make -j8 CFD_solver_gpu"
    echo "  4. Run tests from TESTING_AND_VALIDATION_PLAN.md"
    exit 0
elif [ $FAIL_COUNT -le 2 ] && [ $WARN_COUNT -eq 0 ]; then
    echo -e "${YELLOW}⚠ MINOR ISSUES FOUND${NC}"
    echo "Review logs in: $BUILD_DIR/"
    exit 1
else
    echo -e "${RED}✗ VALIDATION FAILED${NC}"
    echo "Fix errors before proceeding."
    echo "Review logs in: $BUILD_DIR/"
    exit 2
fi

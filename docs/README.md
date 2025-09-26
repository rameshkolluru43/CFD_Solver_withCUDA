# Documentation

This directory contains documentation for the CFD Solver with CUDA project.

## Structure

- **html/** - Generated HTML documentation (created by Doxygen)
- **theory.md** - Mathematical background and theoretical foundation
- **images/** - Images and diagrams used in documentation

## Generating Documentation

### Using the Build Script
```bash
# Generate HTML documentation
./scripts/build_docs.sh

# Generate HTML and PDF documentation
./scripts/build_docs.sh --pdf

# Generate and open documentation in browser
./scripts/build_docs.sh --open
```

### Using Doxygen Directly
```bash
# Generate documentation
doxygen Doxyfile

# View generated documentation
open docs/html/index.html
```

### Using Make/CMake
```bash
# With Make
make docs

# With CMake
cd build
make docs
```

## Documentation Standards

### Code Documentation
All public APIs must be documented using Doxygen format:

```cpp
/**
 * @brief Brief description of the function
 * @param param1 Description of parameter 1
 * @param param2 Description of parameter 2
 * @return Description of return value
 * @throws std::exception When exception is thrown
 * 
 * Detailed description of the function behavior.
 * 
 * @example
 * ```cpp
 * MyClass obj;
 * bool result = obj.myFunction(param1, param2);
 * ```
 */
bool myFunction(int param1, double param2);
```

### CUDA Kernel Documentation
CUDA kernels require additional information:

```cpp
/**
 * @brief CUDA kernel for specific computation
 * @param input_d Device pointer to input data [size: N]
 * @param output_d Device pointer to output data [size: N]
 * @param N Number of elements to process
 * 
 * Launch configuration: Use 1D blocks of 256 threads
 * Memory requirements: 2*N*sizeof(double) device memory
 * Performance: ~100 GB/s on RTX 3080
 */
__global__ void myKernel(const double* input_d, double* output_d, int N);
```

### File Documentation
Each source file should have a header comment:

```cpp
/**
 * @file filename.cpp
 * @brief Brief description of file contents
 * @author Author Name
 * @date Creation date
 * @version Version number
 * 
 * Detailed description of the file's purpose and contents.
 */
```

## Documentation Tools

### Required Tools
- **Doxygen** (1.8.0+) - Documentation generator
- **Graphviz** (optional) - For generating diagrams

### Installation
```bash
# Ubuntu/Debian
sudo apt install doxygen graphviz

# macOS
brew install doxygen graphviz

# CentOS/RHEL
sudo yum install doxygen graphviz
```

## Documentation Deployment

Documentation is automatically built and deployed using GitHub Actions:

1. **Push to main branch** - Triggers documentation build
2. **GitHub Actions** - Runs Doxygen and generates HTML
3. **GitHub Pages** - Hosts the documentation
4. **Access** - Available at https://rameshkolluru43.github.io/CFD_Solver_withCUDA/

## Troubleshooting

### Common Issues

**Doxygen warnings about obsolete tags:**
- Update Doxyfile using: `doxygen -u Doxyfile`

**Missing CUDA include paths:**
- Install CUDA Toolkit or update INCLUDE_PATH in Doxyfile

**Graphviz not found:**
- Install Graphviz package for diagram generation

**LaTeX errors (PDF generation):**
- Install LaTeX distribution (texlive-full on Ubuntu)

### Manual Fixes

If documentation generation fails:

1. Check Doxygen version: `doxygen --version`
2. Verify Doxyfile syntax: `doxygen -g test_config && diff Doxyfile test_config`
3. Check for missing dependencies
4. Review Doxygen output for specific errors

## Contributing to Documentation

1. **Update code comments** when modifying functions
2. **Add examples** for new features
3. **Update theory.md** for new algorithms
4. **Test documentation** generation before submitting PR
5. **Follow style guidelines** for consistency

See [CONTRIBUTING.md](../CONTRIBUTING.md) for detailed contribution guidelines.
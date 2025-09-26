# Contributing to CFD Solver with CUDA

Thank you for your interest in contributing to the CFD Solver with CUDA project! This document provides guidelines and information for contributors.

## 🤝 Ways to Contribute

- **Bug Reports**: Submit detailed bug reports with reproducible examples
- **Feature Requests**: Propose new features or enhancements
- **Code Contributions**: Submit pull requests for bug fixes or new features
- **Documentation**: Improve documentation, tutorials, or examples
- **Testing**: Help expand test coverage or improve existing tests
- **Performance Optimization**: Optimize CUDA kernels or algorithms

## 📋 Before You Start

1. **Check existing issues** to see if your bug/feature has already been reported
2. **Read the documentation** to understand the project structure and goals
3. **Set up your development environment** following the installation guide
4. **Run existing tests** to ensure your environment is working correctly

## 🛠️ Development Setup

### Prerequisites
- CUDA Toolkit 11.0+
- CMake 3.18+
- C++14 compatible compiler
- Git
- Doxygen (for documentation)

### Getting Started
```bash
# Fork the repository on GitHub
# Clone your fork
git clone https://github.com/YOUR_USERNAME/CFD_Solver_withCUDA.git
cd CFD_Solver_withCUDA

# Add upstream remote
git remote add upstream https://github.com/rameshkolluru43/CFD_Solver_withCUDA.git

# Create development branch
git checkout -b feature/your-feature-name

# Build the project
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# Run tests
make test
```

## 📝 Coding Standards

### Code Style
- Follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)
- Use consistent indentation (4 spaces, no tabs)
- Maximum line length: 100 characters
- Use descriptive variable and function names

### Naming Conventions
- **Classes**: PascalCase (`CFDSolver`, `FluidProperties`)
- **Functions**: camelCase (`computeDivergence`, `applyBoundaryConditions`)
- **Variables**: snake_case (`grid_spacing`, `time_step`)
- **Constants**: UPPER_SNAKE_CASE (`MAX_ITERATIONS`, `DEFAULT_TOLERANCE`)
- **Files**: snake_case (`cfd_solver.h`, `cuda_kernels.cu`)

### Documentation
All public APIs must be documented using Doxygen format:

```cpp
/**
 * @brief Brief description of the function
 * @param param1 Description of parameter 1
 * @param param2 Description of parameter 2
 * @return Description of return value
 * @throws std::exception Description of when exception is thrown
 * 
 * Detailed description of the function behavior, including
 * any important implementation details or usage notes.
 * 
 * @example
 * ```cpp
 * CFDSolver solver;
 * solver.initialize(params);
 * solver.run();
 * ```
 */
bool function_name(int param1, double param2);
```

### CUDA Specific Guidelines
- Always check CUDA errors using the `CUDA_CHECK` macro
- Use appropriate memory access patterns for coalescing
- Document kernel launch configurations
- Include performance considerations in comments

```cpp
/**
 * @brief CUDA kernel for computing divergence
 * @param u_d Device pointer to x-velocity [size: nx*ny*nz]
 * @param div_d Device pointer to output divergence [size: nx*ny*nz]
 * @param nx,ny,nz Grid dimensions
 * @param dx,dy,dz Grid spacing
 * 
 * Kernel uses 2D thread blocks of size 16x16 for optimal memory coalescing.
 * Each thread computes divergence for one grid point using central differences.
 * 
 * Performance: ~100 GB/s memory bandwidth on RTX 3080
 */
__global__ void computeDivergence(const double* u_d, double* div_d, ...);
```

## 🧪 Testing

### Test Requirements
- All new features must include unit tests
- Bug fixes must include regression tests
- Tests should cover edge cases and error conditions
- Maintain >80% code coverage

### Test Structure
```cpp
#include <gtest/gtest.h>
#include "cfd_solver.h"

class CFDSolverTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize test data
    }
    
    void TearDown() override {
        // Cleanup
    }
    
    cfd::CFDSolver solver;
    cfd::SimulationParameters params;
};

TEST_F(CFDSolverTest, InitializationSuccess) {
    EXPECT_TRUE(solver.initialize(params));
    EXPECT_GT(solver.getMemoryUsage(), 0);
}
```

### Running Tests
```bash
# Run all tests
make test

# Run specific test suite
./bin/unit_tests --gtest_filter="CFDSolverTest.*"

# Run with coverage
make coverage
```

## 📊 Performance Guidelines

### Optimization Priorities
1. **Correctness**: Ensure numerical accuracy first
2. **Memory bandwidth**: Optimize memory access patterns
3. **Compute utilization**: Maximize GPU occupancy
4. **Host-device transfers**: Minimize data movement

### Benchmarking
- Include performance benchmarks for new kernels
- Compare against CPU baseline when applicable
- Document performance characteristics in code comments

## 🔄 Pull Request Process

### Before Submitting
- [ ] Code follows style guidelines
- [ ] All tests pass
- [ ] Documentation is updated
- [ ] Performance impact is considered
- [ ] Commit messages are descriptive

### PR Template
```
## Description
Brief description of changes and motivation.

## Type of Change
- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update

## Testing
- [ ] Unit tests added/updated
- [ ] Integration tests pass
- [ ] Performance benchmarks included (if applicable)

## Checklist
- [ ] Code follows project style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No new warnings introduced
```

### Review Process
1. **Automated checks**: CI/CD pipeline runs tests and checks
2. **Code review**: Maintainers review code for quality and correctness
3. **Performance review**: Performance impact is evaluated
4. **Documentation review**: Documentation completeness is checked
5. **Final approval**: Maintainer approves and merges

## 🐛 Bug Reports

### Bug Report Template
```
**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Environment:**
- OS: [e.g. Ubuntu 20.04]
- CUDA Version: [e.g. 11.8]
- GPU: [e.g. RTX 3080]
- Compiler: [e.g. GCC 9.4]

**Additional context**
Add any other context about the problem here.
```

## 💡 Feature Requests

### Feature Request Template
```
**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is.

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions.

**Additional context**
Add any other context or screenshots about the feature request here.
```

## 📚 Documentation Guidelines

### Types of Documentation
- **API Documentation**: Doxygen comments in header files
- **User Guide**: README.md and usage examples
- **Theory Documentation**: Mathematical background in docs/theory.md
- **Developer Guide**: Implementation details and architecture

### Documentation Standards
- Use clear, concise language
- Include code examples where appropriate
- Keep documentation up-to-date with code changes
- Use proper markdown formatting

## 🏷️ Release Process

### Version Numbering
We follow semantic versioning (MAJOR.MINOR.PATCH):
- **MAJOR**: Incompatible API changes
- **MINOR**: Backward-compatible functionality additions
- **PATCH**: Backward-compatible bug fixes

### Release Checklist
- [ ] All tests pass
- [ ] Documentation updated
- [ ] Version numbers updated
- [ ] Changelog updated
- [ ] Performance benchmarks run
- [ ] Release notes prepared

## 📞 Getting Help

- **Discussions**: Use GitHub Discussions for questions and ideas
- **Issues**: Use GitHub Issues for bug reports and feature requests
- **Email**: Contact maintainers directly for sensitive issues
- **Documentation**: Check existing documentation first

## 📄 License

By contributing to this project, you agree that your contributions will be licensed under the GNU General Public License v3.0.

## 🙏 Recognition

Contributors are recognized in:
- Project README.md
- Release notes
- Documentation acknowledgments
- Git commit history

Thank you for contributing to CFD Solver with CUDA! 🚀
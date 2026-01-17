# Makefile for Incompressible Flow Solver
# Uses cell-centered staggered grid with finite volume method

# Compiler settings
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra
DEBUG_FLAGS = -g -DDEBUG -O0
INCLUDES = -I../include -I../Basic_Function_Files -I../src

# Library settings - adjust for macOS
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS - use Homebrew paths if available
    LIBS = -ljsoncpp -lm
    HOMEBREW_PREFIX = $(shell brew --prefix 2>/dev/null || echo /usr/local)
    INCLUDES += -I$(HOMEBREW_PREFIX)/include
    LIBS += -L$(HOMEBREW_PREFIX)/lib
else
    # Linux
    LIBS = -ljsoncpp -lm -fopenmp
    CXXFLAGS += -fopenmp
endif

# CUDA settings (optional)
NVCC = nvcc
CUDA_FLAGS = -O3 -std=c++11 -arch=sm_60
CUDA_INCLUDES = -I../CUDA_KERNELS
CUDA_LIBS = -lcudart

# Source directories
SRC_DIR = .
GRID_SRC_DIR = ../src
BASIC_SRC_DIR = ../Basic_Function_Files

# Source files
SOLVER_SOURCES = Incompressible_Solver.cpp \
                Incompressible_Utilities.cpp

MAIN_SOURCES = Incompressible_Main.cpp
TEST_SOURCES = test_standalone.cpp

# Temporarily exclude external dependencies for standalone build
GRID_SOURCES = 
BASIC_SOURCES = 

# Object files
SOLVER_OBJECTS = $(SOLVER_SOURCES:.cpp=.o)
GRID_OBJECTS = $(notdir $(GRID_SOURCES:.cpp=.o))
BASIC_OBJECTS = $(notdir $(BASIC_SOURCES:.cpp=.o))

ALL_OBJECTS = $(SOLVER_OBJECTS) $(GRID_OBJECTS) $(BASIC_OBJECTS)

# Executable names
EXECUTABLE = incompressible_solver
DEBUG_EXECUTABLE = incompressible_solver_debug
TEST_EXECUTABLE = test_standalone
OLD_TARGET = ns_solver

# Default target
all: $(TEST_EXECUTABLE)

# Main target (for full integration)
main: $(EXECUTABLE)

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(DEBUG_EXECUTABLE)

# Object files for standalone test
standalone_test.o: test_standalone.cpp
	$(CXX) $(CXXFLAGS) -c test_standalone.cpp -o standalone_test.o

incompressible_globals.o: incompressible_globals.cpp
	$(CXX) $(CXXFLAGS) -c incompressible_globals.cpp -o incompressible_globals.o

incompressible_solver_standalone.o: Incompressible_Solver_Standalone.cpp
	$(CXX) $(CXXFLAGS) -c Incompressible_Solver_Standalone.cpp -o incompressible_solver_standalone.o

test_standalone: standalone_test.o incompressible_globals.o incompressible_solver_standalone.o
	$(CXX) $(CXXFLAGS) standalone_test.o incompressible_globals.o incompressible_solver_standalone.o -o test_standalone

# Main executable (requires external dependencies)
$(EXECUTABLE): $(ALL_OBJECTS) Incompressible_Main.o
	@echo "Linking $(EXECUTABLE)..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo "Build complete: $(EXECUTABLE)"

# Debug executable
$(DEBUG_EXECUTABLE): $(ALL_OBJECTS)
	@echo "Linking debug executable..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	@echo "Debug build complete: $(DEBUG_EXECUTABLE)"

# Compile solver source files
%.o: %.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile grid source files
$(GRID_OBJECTS): $(GRID_SOURCES)
	@echo "Compiling grid sources..."
	@for src in $(GRID_SOURCES); do \
		if [ -f "$$src" ]; then \
			obj=$$(basename $$src .cpp).o; \
			echo "  Compiling $$src -> $$obj"; \
			$(CXX) $(CXXFLAGS) $(INCLUDES) -c $$src -o $$obj; \
		fi; \
	done

# Compile basic function files  
$(BASIC_OBJECTS): $(BASIC_SOURCES)
	@echo "Compiling basic function sources..."
	@for src in $(BASIC_SOURCES); do \
		if [ -f "$$src" ]; then \
			obj=$$(basename $$src .cpp).o; \
			echo "  Compiling $$src -> $$obj"; \
			$(CXX) $(CXXFLAGS) $(INCLUDES) -c $$src -o $$obj; \
		fi; \
	done

# Legacy target for backward compatibility
$(OLD_TARGET): $(EXECUTABLE)
	ln -sf $(EXECUTABLE) $(OLD_TARGET)

# Clean targets
clean:
	@echo "Cleaning object files..."
	rm -f *.o
	@echo "Clean complete."

distclean: clean
	@echo "Cleaning executables..."
	rm -f $(EXECUTABLE) $(DEBUG_EXECUTABLE) $(TEST_EXECUTABLE) $(OLD_TARGET)
	rm -f test_incompressible_solver.cpp
	rm -rf output/
	@echo "Distclean complete."

# Create output directory
output:
	@mkdir -p output
	@echo "Output directory created."

# Run validation test
validate: $(EXECUTABLE) output
	@echo "Running validation test case..."
	./$(EXECUTABLE) -v > output/validation.log 2>&1
	@echo "Validation complete. Check output/validation.log"

# Run example with lid-driven cavity
example-cavity: $(EXECUTABLE) output
	@echo "Running lid-driven cavity example..."
	./$(EXECUTABLE) -v -o output/cavity/ > output/cavity.log 2>&1
	@echo "Cavity example complete. Check output/cavity/"

# Test run target
test: $(TEST_EXECUTABLE)
	@echo "Running standalone test..."
	./$(TEST_EXECUTABLE)

# Legacy run target
run: $(TEST_EXECUTABLE)
	@echo "Running incompressible solver test..."
	./$(TEST_EXECUTABLE)

# Show help
help:
	@echo "Available targets:"
	@echo "  all              - Build release version (default)"
	@echo "  debug            - Build debug version"
	@echo "  clean            - Remove object files"
	@echo "  distclean        - Remove all generated files"
	@echo "  validate         - Run validation test case"
	@echo "  example-cavity   - Run lid-driven cavity example"
	@echo "  run              - Run solver with default settings"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Usage examples:"
	@echo "  make                    # Build release version"
	@echo "  make debug              # Build debug version"
	@echo "  make validate           # Run validation"
	@echo "  make clean && make all  # Clean rebuild"

# Phony targets
.PHONY: all debug clean distclean output validate example-cavity run help
# Compiler and flags
CXX = mpic++
CXXFLAGS = -O2 -Wall
TARGET = ns_solver
SRC = Incompressible_NS.cpp

# Build target
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean up build files
clean:
	rm -f $(TARGET) *.vtk

# Run with specified number of processes
run: $(TARGET)
	mpirun -np 4 ./$(TARGET)

.PHONY: all clean run
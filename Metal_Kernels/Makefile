# Apple Silicon / macOS only: Objective-C++ bridge + Metal
SDKROOT ?= $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)
CXX = clang++
CXXFLAGS = -std=c++17 -O2 -fobjc-arc -Iinclude \
	-isysroot $(SDKROOT) \
	-F$(SDKROOT)/System/Library/Frameworks
LDFLAGS = -framework Metal -framework Foundation

SRC = src/MetalHostBridge.mm
OBJ = build/MetalHostBridge.o
LIB = build/libmetal_cfd_bridge.a

.PHONY: all clean bridge metallib

all: bridge

bridge: $(LIB)

build:
	mkdir -p build

$(OBJ): $(SRC) include/MetalKernelBridge.h | build
	$(CXX) $(CXXFLAGS) -c $(SRC) -o $(OBJ)

$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)

metallib:
	./scripts/build_metallib.sh

clean:
	rm -rf build compiled

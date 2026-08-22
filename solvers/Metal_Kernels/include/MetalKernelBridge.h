#pragma once

#include <cstddef>
#include <cstdint>

struct MetalCFDContext;

// Owns Metal device, loaded metallib, cached pipeline states, command queue.
MetalCFDContext* metal_cfd_create(const char* metallib_path, char* err_buf, size_t err_len);
void metal_cfd_destroy(MetalCFDContext* ctx);

// Returns id<MTLComputePipelineState> as opaque pointer; nil on failure.
void* metal_cfd_pipeline(MetalCFDContext* ctx, const char* kernel_name);

// Function names in CFDSolverKernels.metal (use with newFunctionWithName:).
extern const char kMetalKernelEulerExplicitUpdate4[];
extern const char kMetalKernelEulerExplicitUpdateInplace4[];
extern const char kMetalKernelRusanovFaceFlux2D[];
extern const char kMetalKernelInviscidDtDenominatorQuad[];
extern const char kMetalKernelConservativeToPrimitive2D[];

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include "MetalKernelBridge.h"

#include <cstring>

const char kMetalKernelEulerExplicitUpdate4[] = "cfd_euler_explicit_update_4";
const char kMetalKernelEulerExplicitUpdateInplace4[] = "cfd_euler_explicit_update_inplace_4";
const char kMetalKernelRusanovFaceFlux2D[] = "cfd_rusanov_face_flux_2d";
const char kMetalKernelInviscidDtDenominatorQuad[] = "cfd_inviscid_dt_denominator_quad";
const char kMetalKernelConservativeToPrimitive2D[] = "cfd_conservative_to_primitive_2d";

@interface MetalCFDContextInternal : NSObject
@property(nonatomic, strong) id<MTLDevice> device;
@property(nonatomic, strong) id<MTLLibrary> library;
@property(nonatomic, strong) id<MTLCommandQueue> queue;
@property(nonatomic, strong) NSMutableDictionary<NSString*, id<MTLComputePipelineState>>* pipelines;
@end

@implementation MetalCFDContextInternal
@end

static void copy_err(char* err_buf, size_t err_len, NSString* msg)
{
    if (!err_buf || err_len == 0)
        return;
    const char* c = msg.UTF8String;
    std::strncpy(err_buf, c ? c : "unknown error", err_len - 1);
    err_buf[err_len - 1] = '\0';
}

MetalCFDContext* metal_cfd_create(const char* metallib_path, char* err_buf, size_t err_len)
{
    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev) {
        copy_err(err_buf, err_len, @"No Metal device");
        return nullptr;
    }
    NSError* err = nil;
    NSURL* url = [NSURL fileURLWithPath:[NSString stringWithUTF8String:metallib_path]];
    id<MTLLibrary> lib = [dev newLibraryWithURL:url error:&err];
    if (!lib) {
        copy_err(err_buf, err_len, err.localizedDescription ?: @"Failed to load metallib");
        return nullptr;
    }
    MetalCFDContextInternal* box = [MetalCFDContextInternal new];
    box.device = dev;
    box.library = lib;
    box.queue = [dev newCommandQueue];
    box.pipelines = [NSMutableDictionary dictionary];
    return (MetalCFDContext*)CFBridgingRetain(box);
}

void metal_cfd_destroy(MetalCFDContext* ctx)
{
    if (!ctx)
        return;
    MetalCFDContextInternal* box
        = (__bridge_transfer MetalCFDContextInternal*)((void*)ctx);
    (void)box;
}

void* metal_cfd_pipeline(MetalCFDContext* ctx, const char* kernel_name)
{
    if (!ctx || !kernel_name)
        return nullptr;
    MetalCFDContextInternal* box = (__bridge MetalCFDContextInternal*)((void*)ctx);
    NSString* key = [NSString stringWithUTF8String:kernel_name];
    id<MTLComputePipelineState> ps = box.pipelines[key];
    if (ps)
        return (__bridge void*)ps;
    NSError* err = nil;
    id<MTLFunction> fn = [box.library newFunctionWithName:key];
    if (!fn)
        return nullptr;
    ps = [box.device newComputePipelineStateWithFunction:fn error:&err];
    if (!ps)
        return nullptr;
    box.pipelines[key] = ps;
    return (__bridge void*)ps;
}

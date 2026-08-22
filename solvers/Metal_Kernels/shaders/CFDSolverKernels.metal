//
// CFDSolverKernels.metal
// Metal compute kernels analogous to CUDA_KERNELS/Time_Integration_Cuda_Kernels.cu
// and Rusanov face flux (2D Euler, four conservative variables).
//
// Build: ../scripts/build_metallib.sh
//

#include <metal_stdlib>
using namespace metal;

static constexpr float kGamma = 1.4f;
static constexpr float kGammaM1 = kGamma - 1.0f;

inline void cfd_primitive_from_U(float4 U, thread float& rho, thread float& u, thread float& v,
    thread float& P, thread float& c)
{
    rho = U.x;
    float invR = 1.0f / max(rho, 1e-20f);
    u = U.y * invR;
    v = U.z * invR;
    float ke = 0.5f * (u * u + v * v);
    P = kGammaM1 * (U.w - rho * ke);
    P = max(P, 1e-20f);
    c = sqrt(kGamma * P * invR);
}

kernel void cfd_euler_explicit_update_4(device const float* U_old [[buffer(0)]],
    device const float* net_flux [[buffer(1)]], device float* U_new [[buffer(2)]],
    device const float* inv_area [[buffer(3)]], constant float& dt [[buffer(4)]],
    constant uint& n_cells [[buffer(5)]], uint gid [[thread_position_in_grid]])
{
    if (gid >= n_cells)
        return;
    size_t base = static_cast<size_t>(gid) * 4;
    float ia = inv_area[gid];
    for (int k = 0; k < 4; k++)
        U_new[base + k] = U_old[base + k] - dt * ia * net_flux[base + k];
}

kernel void cfd_euler_explicit_update_inplace_4(device float* U [[buffer(0)]],
    device const float* net_flux [[buffer(1)]], device const float* inv_area [[buffer(2)]],
    constant float& dt [[buffer(3)]], constant uint& n_cells [[buffer(4)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= n_cells)
        return;
    size_t base = static_cast<size_t>(gid) * 4;
    float ia = inv_area[gid];
    for (int k = 0; k < 4; k++)
        U[base + k] -= dt * ia * net_flux[base + k];
}

kernel void cfd_rusanov_face_flux_2d(device const float* U_L [[buffer(0)]],
    device const float* U_R [[buffer(1)]], device const float* face_n [[buffer(2)]],
    device const float* face_dl [[buffer(3)]], device float* flux_out [[buffer(4)]],
    constant uint& n_faces [[buffer(5)]], uint gid [[thread_position_in_grid]])
{
    if (gid >= n_faces)
        return;
    size_t b = static_cast<size_t>(gid) * 4;
    float4 UL(U_L[b], U_L[b + 1], U_L[b + 2], U_L[b + 3]);
    float4 UR(U_R[b], U_R[b + 1], U_R[b + 2], U_R[b + 3]);

    float rL, uL, vL, pL, cL;
    float rR, uR, vR, pR, cR;
    cfd_primitive_from_U(UL, rL, uL, vL, pL, cL);
    cfd_primitive_from_U(UR, rR, uR, vR, pR, cR);

    float nx = face_n[gid * 2];
    float ny = face_n[gid * 2 + 1];
    float dl = face_dl[gid];

    float vnL = uL * nx + vL * ny;
    float vnR = uR * nx + vR * ny;

    float FL0 = rL * vnL * dl;
    float FL1 = (rL * uL * vnL + pL * nx) * dl;
    float FL2 = (rL * vL * vnL + pL * ny) * dl;
    float FL3 = vnL * (UL.w + pL) * dl;

    float FR0 = rR * vnR * dl;
    float FR1 = (rR * uR * vnR + pR * nx) * dl;
    float FR2 = (rR * vR * vnR + pR * ny) * dl;
    float FR3 = vnR * (UR.w + pR) * dl;

    float lam = max(fabs(vnL) + cL, fabs(vnR) + cR);
    float lamdl = lam * dl;

    flux_out[b] = 0.5f * (FL0 + FR0) - 0.5f * lamdl * (UR.x - UL.x);
    flux_out[b + 1] = 0.5f * (FL1 + FR1) - 0.5f * lamdl * (U_R[b + 1] - U_L[b + 1]);
    flux_out[b + 2] = 0.5f * (FL2 + FR2) - 0.5f * lamdl * (U_R[b + 2] - U_L[b + 2]);
    flux_out[b + 3] = 0.5f * (FL3 + FR3) - 0.5f * lamdl * (U_R[b + 3] - U_L[b + 3]);
}

kernel void cfd_inviscid_dt_denominator_quad(device const float* face_pack [[buffer(0)]],
    device const float* rho [[buffer(1)]], device const float* u [[buffer(2)]],
    device const float* v [[buffer(3)]], device const float* P [[buffer(4)]],
    device float* denom_out [[buffer(5)]], constant uint& n_cells [[buffer(6)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= n_cells)
        return;
    float r = max(rho[gid], 1e-20f);
    float ui = u[gid];
    float vi = v[gid];
    float pi = max(P[gid], 1e-20f);
    float c = sqrt(kGamma * pi / r);
    size_t o = static_cast<size_t>(gid) * 12;
    float s = 0.0f;
    for (int f = 0; f < 4; f++)
    {
        float nx = face_pack[o + f * 3];
        float ny = face_pack[o + f * 3 + 1];
        float dl = face_pack[o + f * 3 + 2];
        float vn = ui * nx + vi * ny;
        s += (fabs(vn) + c) * dl;
    }
    denom_out[gid] = s;
}

kernel void cfd_conservative_to_primitive_2d(device const float* U [[buffer(0)]],
    device float* rho [[buffer(1)]], device float* u [[buffer(2)]], device float* v [[buffer(3)]],
    device float* P [[buffer(4)]], constant uint& n_cells [[buffer(5)]], uint gid [[thread_position_in_grid]])
{
    if (gid >= n_cells)
        return;
    size_t b = static_cast<size_t>(gid) * 4;
    float4 Uv(U[b], U[b + 1], U[b + 2], U[b + 3);
    float r, ux, vy, p, cc;
    cfd_primitive_from_U(Uv, r, ux, vy, p, cc);
    rho[gid] = r;
    u[gid] = ux;
    v[gid] = vy;
    P[gid] = p;
}

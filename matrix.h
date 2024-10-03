/*
 SIMD++/matrix.h
 
 Created by Treata Norouzi on 9/20/24.
 
 See LICENSE for the licensing information.
*/

#include <simd/matrix.h>

#define SIMD_DEFINE_VEC_TYPE(type, n) \
\
float type##n##_len(type##n const v) \
{ \
    return sqrtf(simd_dot(v, v)); \
} \
\
void type##n##_min(type##n r, type##n const a, type##n const b) \
{ \
    int i; \
    for (i = 0; i < n; ++i) \
        r[i] = a[i] < b[i] ? a[i] : b[i]; \
} \
\
void type##n##_max(type##n r, type##n const a, type##n const b) \
{ \
    int i; \
    for (i = 0; i < n; ++i) \
        r[i] = a[i] > b[i] ? a[i] : b[i]; \
} \
//\
//void type##n##_dup(type##n r, type##n const src) \
//{ \
//    int i; \
//    for(i=0; i<n; ++i) \
//        r[i] = src[i]; \
//}

#ifdef __cplusplus
extern "C" {
#endif

// - Define SIMD types and their corresponding functions

//SIMD_DEFINE_VEC_TYPE(simd_int, 2)
//SIMD_DEFINE_VEC_TYPE(simd_int, 3)
//SIMD_DEFINE_VEC_TYPE(simd_int, 4)

SIMD_DEFINE_VEC_TYPE(simd_half, 2)
SIMD_DEFINE_VEC_TYPE(simd_half, 3)
SIMD_DEFINE_VEC_TYPE(simd_half, 4)

SIMD_DEFINE_VEC_TYPE(simd_float, 2)
SIMD_DEFINE_VEC_TYPE(simd_float, 3)
SIMD_DEFINE_VEC_TYPE(simd_float, 4)

SIMD_DEFINE_VEC_TYPE(simd_double, 2)
SIMD_DEFINE_VEC_TYPE(simd_double, 3)
SIMD_DEFINE_VEC_TYPE(simd_double, 4)

// MARK: - Vectors

// TODO: Half and Double
//void vec3_mul_cross(simd_float3 r, simd_float3 const a, simd_float3 const b)
//{
//    r[0] = a[1]*b[2] - a[2]*b[1];
//    r[1] = a[2]*b[0] - a[0]*b[2];
//    r[2] = a[0]*b[1] - a[1]*b[0];
//}
//
//void vec4_mul_cross(simd_float4 r, simd_float4 const a, simd_float4 const b)
//{
//    r[0] = a[1]*b[2] - a[2]*b[1];
//    r[1] = a[2]*b[0] - a[0]*b[2];
//    r[2] = a[0]*b[1] - a[1]*b[0];
//    r[3] = 1.f;
//}

// MARK: - Matrices

// TODO: Rename to simd_floatNxN_scale_aniso
/// Anisotropic scaling
simd_half4x4 SIMD_CFUNC simd_scale_aniso(simd_half4x4 M, simd_half4 v);
simd_float4x4 SIMD_CFUNC simd_scale_aniso(simd_float4x4 M, simd_float4 v);
simd_double4x4 SIMD_CFUNC simd_scale_aniso(simd_double4x4 M, simd_double4 v);

simd_half4x4 SIMD_CFUNC simd_translate(simd_half3 v);
simd_float4x4 SIMD_CFUNC simd_translate(simd_float3 v);
simd_double4x4 SIMD_CFUNC simd_translate(simd_double3 v);

// TODO: float4x4_translate_in_place

// TODO: float4x4_from_vec3_mul_outer

// TODO: float4x4_rotate

simd_half4x4 SIMD_CFUNC simd_xRotate(simd_half1 angleRadians);
simd_float4x4 SIMD_CFUNC simd_xRotate(float angleRadians);
simd_double4x4 SIMD_CFUNC simd_xRotate(double angleRadians);

simd_half4x4 SIMD_CFUNC simd_yRotate(simd_half1 angleRadians);
simd_float4x4 SIMD_CFUNC simd_yRotate(float angleRadians);
simd_double4x4 SIMD_CFUNC simd_yRotate(double angleRadians);

simd_half4x4 SIMD_CFUNC simd_zRotate(simd_half1 angleRadians);
simd_float4x4 SIMD_CFUNC simd_zRotate(float angleRadians);
simd_double4x4 SIMD_CFUNC simd_zRotate(double angleRadians);

void SIMD_CFUNC simd_float4x4_orthonormalize(simd_float4x4 R, simd_float4x4 const M);

// TODO: float4x4_frustum

// TODO: float4x4_ortho

// TODO: float4x4_perspective

// TODO: float4x4_look_at

// TODO: quat_conj - already exists?

// TODO: quat_rotate

// TODO: quat_mul_vec3 - already exists?

// TODO: float4x4_from_quat - already exists?

// TODO: float4x4o_mul_quat - already exists?

// TODO: quat_from_float4x4

// TODO: float4x4_arcball

#ifdef __cplusplus
} /* extern "C" */

// MARK: - C++

namespace simd {

    SIMD_CPPFUNC half4x4 scale_aniso(half4x4 M, half4 v) { return ::simd_scale_aniso(M, v); }
    SIMD_CPPFUNC float4x4 scale_aniso(float4x4 M, float v) { return ::simd_scale_aniso(M, v); }
    SIMD_CPPFUNC double4x4 scale_aniso(double4x4 M, double4 v) { return ::simd_scale_aniso(M, v); }

    SIMD_CPPFUNC half4x4 translate(const half3 &v) { return ::simd_translate(v); }
    SIMD_CPPFUNC float4x4 translate(const float3 &v) { return ::simd_translate(v); }
    SIMD_CPPFUNC double4x4 translate(const double3 &v) { return ::simd_translate(v); }

    // SIMD_CPPFUNC float4x4 rotate(float4x4 M, const float3 v, float angle) { return ::simd_rotate(M, v, angle); }

    SIMD_CPPFUNC half4x4 xRotate(simd_half1 angleRadians) { return ::simd_xRotate(angleRadians); }
    SIMD_CPPFUNC float4x4 xRotate(float angleRadians) { return ::simd_xRotate(angleRadians); }
    SIMD_CPPFUNC double4x4 xRotate(double angleRadians) { return ::simd_xRotate(angleRadians); }

    SIMD_CPPFUNC half4x4 yRotate(simd_half1 angleRadians) { return ::simd_yRotate(angleRadians); }
    SIMD_CPPFUNC float4x4 yRotate(float angleRadians) { return ::simd_yRotate(angleRadians); }
    SIMD_CPPFUNC double4x4 yRotate(double angleRadians) { return ::simd_yRotate(angleRadians); }

    SIMD_CPPFUNC half4x4 zRotate(simd_half1 angleRadians) { return ::simd_zRotate(angleRadians); }
    SIMD_CPPFUNC float4x4 zRotate(float angleRadians) { return ::simd_zRotate(angleRadians); }
    SIMD_CPPFUNC double4x4 zRotate(double angleRadians) { return ::simd_zRotate(angleRadians); }

    SIMD_CPPFUNC void simd_float4x4_orthonormalize(simd_float4x4 R, simd_float4x4 const M) { return ::simd_float4x4_orthonormalize(R, M); };

}

extern "C" {
#endif /* __cplusplus */

#pragma mark - C Implementation

simd_half4x4 SIMD_CFUNC simd_scale_aniso(simd_half4x4 M, simd_half4 v)
{
    return simd_mul(M, simd_diagonal_matrix(v));
}
simd_float4x4 SIMD_CFUNC simd_scale_aniso(simd_float4x4 M, simd_float4 v)
{
    return simd_mul(M, simd_diagonal_matrix(v));
}
simd_double4x4 SIMD_CFUNC simd_scale_aniso(simd_double4x4 M, simd_double4 v)
{
    return simd_mul(M, simd_diagonal_matrix(v));
}

simd_half4x4 SIMD_CFUNC simd_translate(simd_half3 v)
{
    const simd_half4 col0 = { 1.0f, 0.0f, 0.0f, 0.0f };
    const simd_half4 col1 = { 0.0f, 1.0f, 0.0f, 0.0f };
    const simd_half4 col2 = { 0.0f, 0.0f, 1.0f, 0.0f };
    const simd_half4 col3 = { v.x,   v.y,  v.z, 1.0f };
    return simd_matrix(col0, col1, col2, col3);
}
simd_float4x4 SIMD_CFUNC simd_translate(simd_float3 v)
{
    const simd_float4 col0 = { 1.0f, 0.0f, 0.0f, 0.0f };
    const simd_float4 col1 = { 0.0f, 1.0f, 0.0f, 0.0f };
    const simd_float4 col2 = { 0.0f, 0.0f, 1.0f, 0.0f };
    const simd_float4 col3 = { v.x,   v.y,  v.z, 1.0f };
    return simd_matrix(col0, col1, col2, col3);
}
simd_double4x4 SIMD_CFUNC simd_translate(simd_double3 v)
{
    const simd_double4 col0 = { 1.0f, 0.0f, 0.0f, 0.0f };
    const simd_double4 col1 = { 0.0f, 1.0f, 0.0f, 0.0f };
    const simd_double4 col2 = { 0.0f, 0.0f, 1.0f, 0.0f };
    const simd_double4 col3 = { v.x,   v.y,  v.z, 1.0f };
    return simd_matrix(col0, col1, col2, col3);
}



simd_half4x4 SIMD_CFUNC simd_xRotate(simd_half1 angleRadians)
{
    const simd_half1 s = sinf(angleRadians);
    const simd_half1 c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_half4){ 1.0f, 0.0f, 0.0f, 0.0f },
                                 (simd_half4){ 0.0f,    c,    s, 0.0f },
                                 (simd_half4){ 0.0f,   -s,    c, 0.0f },
                                 (simd_half4){ 0.0f, 0.0f, 0.0f, 1.0f });
}
simd_float4x4 SIMD_CFUNC simd_xRotate(float angleRadians)
{
    const float s = sinf(angleRadians);
    const float c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_float4){ 1.0f, 0.0f, 0.0f, 0.0f },
                                 (simd_float4){ 0.0f,    c,    s, 0.0f },
                                 (simd_float4){ 0.0f,   -s,    c, 0.0f },
                                 (simd_float4){ 0.0f, 0.0f, 0.0f, 1.0f });
}
simd_double4x4 SIMD_CFUNC simd_xRotate(double angleRadians)
{
    const double s = sin(angleRadians);
    const double c = cos(angleRadians);
    return simd_matrix_from_rows((simd_double4){ 1.0f, 0.0f, 0.0f, 0.0f },
                                 (simd_double4){ 0.0f,    c,    s, 0.0f },
                                 (simd_double4){ 0.0f,   -s,    c, 0.0f },
                                 (simd_double4){ 0.0f, 0.0f, 0.0f, 1.0f });
}

simd_half4x4 SIMD_CFUNC simd_yRotate(simd_half1 angleRadians)
{
    const simd_half1 s = sinf(angleRadians);
    const simd_half1 c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_half4){    c, 0.0f,    s, 0.0f },
                                 (simd_half4){ 0.0f, 1.0f, 0.0f, 0.0f },
                                 (simd_half4){   -s, 0.0f,    c, 0.0f },
                                 (simd_half4){ 0.0f, 0.0f, 0.0f, 1.0f });
}
simd_float4x4 SIMD_CFUNC simd_yRotate(float angleRadians)
{
    const float s = sinf(angleRadians);
    const float c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_float4){    c, 0.0f,    s, 0.0f },
                                 (simd_float4){ 0.0f, 1.0f, 0.0f, 0.0f },
                                 (simd_float4){   -s, 0.0f,    c, 0.0f },
                                 (simd_float4){ 0.0f, 0.0f, 0.0f, 1.0f });
}
simd_double4x4 SIMD_CFUNC simd_yRotate(double angleRadians)
{
    const double s = sinf(angleRadians);
    const double c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_double4){    c, 0.0f,    s, 0.0f },
                                 (simd_double4){ 0.0f, 1.0f, 0.0f, 0.0f },
                                 (simd_double4){   -s, 0.0f,    c, 0.0f },
                                 (simd_double4){ 0.0f, 0.0f, 0.0f, 1.0f });
}

simd_half4x4 SIMD_CFUNC simd_zRotate(simd_half1 angleRadians)
{
    const simd_half1 s = sinf(angleRadians);
    const simd_half1 c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_half4){    c,    s, 0.0f, 0.0f },
                                 (simd_half4){   -s,    c, 0.0f, 0.0f },
                                 (simd_half4){ 0.0f, 0.0f, 1.0f, 0.0f },
                                 (simd_half4){ 0.0f, 0.0f, 0.0f, 1.0f });
}
simd_float4x4 SIMD_CFUNC simd_zRotate(float angleRadians)
{
    const float s = sinf(angleRadians);
    const float c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_float4){    c,    s, 0.0f, 0.0f },
                                 (simd_float4){   -s,    c, 0.0f, 0.0f },
                                 (simd_float4){ 0.0f, 0.0f, 1.0f, 0.0f },
                                 (simd_float4){ 0.0f, 0.0f, 0.0f, 1.0f });
}
simd_double4x4 SIMD_CFUNC simd_zRotate(double angleRadians)
{
    const double s = sinf(angleRadians);
    const double c = cosf(angleRadians);
    return simd_matrix_from_rows((simd_double4){    c,    s, 0.0f, 0.0f },
                                 (simd_double4){   -s,    c, 0.0f, 0.0f },
                                 (simd_double4){ 0.0f, 0.0f, 1.0f, 0.0f },
                                 (simd_double4){ 0.0f, 0.0f, 0.0f, 1.0f });
}

// TODO: Implement
void SIMD_CFUNC simd_float4x4_orthonormalize(simd_float4x4 R, simd_float4x4 const M)
{

}

#ifdef __cplusplus
}
#endif

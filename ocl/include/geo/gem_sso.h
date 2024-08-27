#ifndef __geo_gem_sso__
#define __geo_gem_sso__

/*  Library of helper functions to support Geometric Element Deformation algorithms with Stretch/Shrink Operation. Implemented from 
    https://www.researchgate.net/publication/277921687_Smoothing_Algorithm_for_Planar_and_Surface_Mesh_Based_on_Element_Geometric_Deformation
    adapted to OpenCL. Specialized implementation provided for triangle and quads as avoiding global storage arrays is considerably faster.

Methods:
    - project_on_plane_tris: projects all points of a triangle on the plane defined by normal and centroid.
                             Overload available for plane defined with an equation form.
    - project_on_plane_quad: projects all points of a quad on the plane defined by normal and centroid.
                             Overload available for plane defined with an equation form.
    - project_on_plane_ngon: projects all points of a n-gon on the plane defined by normal and centroid.
                             Overload available for plane defined with an equation form.
    - stretch_tris: project each vertex of a triangle outward in the direction perpendicular from the opposite edge.
                    This operation regularize the shape but causes a size increase.
    - stretch_quad: project each vertex of a quad outward in the direction perpendicular from the line passing through the previous and the next vertex.
                    This operation regularize the shape but causes a size increase.
    - stretch_ngon: project each vertex of an n-gon outward in the direction perpendicular from the line passing through the previous and the next vertex.
                    This operation regularize the shape but causes a size increase.
    - shrink_tris: scale the triangle by a shrink factor `k` maintaining the old centroid.
    - shrink_quad: scale the quad by a shrink factor `k` maintaining the old centroid.
    - shrink_ngon: scale the n-gon by a shrink factor `k` maintaining the old centroid.
    - pt_surface_distance: given a point P and the corresponding Quadric Error Metric matrix as a flat array of 10 entries,
                           returns the square distance between the point and the surface.
    - squared_error_tris: returns the total squared distance error for each point of the triangle.
    - squared_error_quad: returns the total squared distance error for each point of the quad.
    - squared_error_ngon: returns the total squared distance error for each point of the n-gon.
    - ternary_search_tris: returns the optimal value of `k` which minimizes the squared error distance
    - ternary_search_quad: returns the optimal value of `k` which minimizes the squared error distance
    - ternary_search_ngon: returns the optimal value of `k` which minimizes the squared error distance
*/

#include "geo/tracing.h"

#define TERNARY_SEARCH_MAX_ITERATION 100

static inline __attribute__((overloadable)) void project_on_plane_tris(double3* A, double3* B, double3* C,
    const double3 prim_centroid, const double3 prim_normal)
{
    *A = project_on_plane(*A, prim_centroid, prim_normal);
    *B = project_on_plane(*B, prim_centroid, prim_normal);
    *C = project_on_plane(*C, prim_centroid, prim_normal);
}

static inline __attribute__((overloadable)) void project_on_plane_tris(double3* A, double3* B, double3* C,
    const double4 abcd)
{
    *A = project_on_plane(*A, abcd);
    *B = project_on_plane(*B, abcd);
    *C = project_on_plane(*C, abcd);
}

static inline __attribute__((overloadable)) void project_on_plane_quad(double3* A, double3* B, double3* C, double3* D,
    const double3 prim_centroid, const double3 prim_normal)
{
    *A = project_on_plane(*A, prim_centroid, prim_normal); 
    *B = project_on_plane(*B, prim_centroid, prim_normal);
    *C = project_on_plane(*C, prim_centroid, prim_normal);
    *D = project_on_plane(*D, prim_centroid, prim_normal);
}

static inline __attribute__((overloadable)) void project_on_plane_quad(double3* A, double3* B, double3* C, double3* D,
    const double4 abcd)
{
    *A = project_on_plane(*A, abcd); 
    *B = project_on_plane(*B, abcd);
    *C = project_on_plane(*C, abcd);
    *D = project_on_plane(*D, abcd);
}

static inline __attribute__((overloadable)) void project_on_plane_ngon(int* primpoints, double* P, double* P_planar, const int size, const int idx,
    const double3 prim_centroid, const double3 prim_normal)
{
    for(int i=idx; i<idx+size; i++)
    {
        vstore3(project_on_plane(vload3(primpoints[i], P), prim_centroid, prim_normal), i, P_planar);
    }
}

static inline __attribute__((overloadable)) void project_on_plane_ngon(int* primpoints, double* P, double* P_planar, const int size, const int idx,
    const double4 abcd)
{
    for(int i=idx; i<idx+size; i++)
    {
        vstore3(project_on_plane(vload3(primpoints[i], P), abcd), i, P_planar);
    }
}

static inline void stretch_tris(const double3 A, const double3 B, const double3 C,
    double3* A_stretch, double3* B_stretch, double3* C_stretch,
    double3* centroid_src, double3* centroid_stretch, const double3 prim_normal, const float lambda)
{
    *centroid_src = (A + B + C) / 3.0;

    double3 N  = (double3)(0.0);
    double3 edge = (double3)(0.0);

    // extrude A
    edge = normalize(B-C);
    N = (cross(prim_normal, edge));
    *A_stretch = A + N * lambda;

    // extrude B
    edge = normalize(C-A);
    N = (cross(prim_normal, edge));
    *B_stretch = B + N * lambda;

    // extrude C
    edge = normalize(A-B);
    N = (cross(prim_normal, edge));
    *C_stretch = C + N * lambda;

    *centroid_stretch = (*A_stretch + *B_stretch + *C_stretch) / 3.0;
}

static inline void stretch_quad(const double3 A, const double3 B, const double3 C, const double3 D,
    double3* A_stretch, double3* B_stretch, double3* C_stretch, double3* D_stretch,
    double3* centroid_src, double3* centroid_stretch, const double3 prim_normal, const float lambda)
{
    *centroid_src = (A + B + C + D) / 4.0;

    double3 N  = (double3)(0.0);
    double3 edge = (double3)(0.0);

    // extrude A
    edge = normalize(B-D);
    N = (cross(prim_normal, edge));
    *A_stretch = A + N * lambda;

    // extrude B
    edge = normalize(C-A);
    N = (cross(prim_normal, edge));
    *B_stretch = B + N * lambda;

    // extrude C
    edge = normalize(D-B);
    N = (cross(prim_normal, edge));
    *C_stretch = C + N * lambda;

    // extrude D
    edge = normalize(A-C);
    N = (cross(prim_normal, edge));
    *D_stretch = D + N * lambda;

    *centroid_stretch = (*A_stretch + *B_stretch + *C_stretch + *D_stretch) / 4.0;
}

static inline void stretch_ngon(double* P, double* P_stretch, const int size, const int idx,
    double3* centroid_src, double3* centroid_stretch, const double3 prim_normal, const float lambda)
{
    double3 N = (double3)(0.0);
    double3 edge = (double3)(0.0);
    double3 P1 = (double3)(0.0);

    // re-initializing due to previous data inside? 
    *centroid_src = (double3)(0.0);
    *centroid_stretch = (double3)(0.0);

    int pt_curr, pt_prev, pt_next = 0;
    double3 P_curr, P_prev, P_next = (double3)(0.0);
    for(int i=0; i<size; i++)
    {
        // grab previous/current/next point indexes.
        pt_curr = i;
        pt_prev = i-1; if(pt_prev<0)     { pt_prev += size; }
        pt_next = i+1; if(pt_next>=size) { pt_next -= size; }

        P_curr = vload3(pt_curr+idx, P);
        P_prev = vload3(pt_prev+idx, P);
        P_next = vload3(pt_next+idx, P);

        edge = normalize(P_next - P_prev);
        N = cross(prim_normal, edge);
        P1 = P_curr + N * lambda;

        *centroid_src += P_curr;
        *centroid_stretch += P1;
        vstore3(P1, pt_curr+idx, P_stretch);
    }

    *centroid_src /= (double)(size);
    *centroid_stretch /= (double)(size);
}

static inline void shrink_tris(const double3 A, const double3 B, const double3 C,
    double3* A_shrink, double3* B_shrink, double3* C_shrink,
    const double k, const double3 centroid_src, const double3 centroid_stretch)
{
    *A_shrink = centroid_src + (A - centroid_stretch) * k;
    *B_shrink = centroid_src + (B - centroid_stretch) * k;
    *C_shrink = centroid_src + (C - centroid_stretch) * k;
}

static inline void shrink_quad(const double3 A, const double3 B, const double3 C, const double3 D,
    double3* A_shrink, double3* B_shrink, double3* C_shrink, double3* D_shrink,
    const double k, const double3 centroid_src, const double3 centroid_stretch)
{
    *A_shrink = centroid_src + (A - centroid_stretch) * k;
    *B_shrink = centroid_src + (B - centroid_stretch) * k;
    *C_shrink = centroid_src + (C - centroid_stretch) * k;
    *D_shrink = centroid_src + (D - centroid_stretch) * k;
}

static inline void shrink_ngon(double* P_stretch, double* P_shrink, const int size, const int idx,
    const double k, const double3 centroid_src, const double3 centroid_stretch)
{
    double3 N_shrink = (double3)(0.0);
    for(int i=idx; i<idx+size; i++)
    {
        N_shrink = centroid_src + (vload3(i, P_stretch) - centroid_stretch) * k;
        vstore3(N_shrink, i, P_shrink);
    }
}

static inline double pt_surface_distance(const double3 P, double* QEM, const int offset)
{
    /* Given a position P and the quadric error matrix QEM, computes the quadric error */

    return QEM[offset+0]*P.x*P.x + 2.0*QEM[offset+1]*P.x*P.y + 2.0*QEM[offset+2]*P.x*P.z + 2.0*QEM[offset+3]*P.x
                                 +     QEM[offset+4]*P.y*P.y + 2.0*QEM[offset+5]*P.y*P.z + 2.0*QEM[offset+6]*P.y
                                                             +     QEM[offset+7]*P.z*P.z + 2.0*QEM[offset+8]*P.z 
                                                                                         +     QEM[offset+9];
}

static inline double squared_error_tris(const double3 A, const double3 B, const double3 C, double* QEM, const int idx)
{
    return pt_surface_distance(A, QEM, idx   )
         + pt_surface_distance(B, QEM, idx+10)
         + pt_surface_distance(C, QEM, idx+20);

}

static inline double squared_error_quad(const double3 A, const double3 B, const double3 C, const double3 D, double* QEM, const int idx)
{
    return pt_surface_distance(A, QEM, idx   )
         + pt_surface_distance(B, QEM, idx+10)
         + pt_surface_distance(C, QEM, idx+20)
         + pt_surface_distance(D, QEM, idx+30);
}

static inline double squared_error_ngon(double* P, double* QEM, const int p_idx, const int q_idx, const int size)
{
    double error = 0.0;
    double3 P_offset = (double3)(0.0);
    for(int i=0; i<size; i++)
    {
        P_offset = vload3(i+p_idx, P);
        error += pt_surface_distance(P_offset, QEM, q_idx+i*10);
    }
    return error;
}

static inline double ternary_search_tris(double min, double max, double precision,
    const double3 A, const double3 B, const double3 C, // stretch triangle
    double* QEM, const int idx, const double3 centroid_src, const double3 centroid_stretch)
{
    int COUNTER = 0;
    double min_third, max_third = 0.0;
    double3 A1, A2, B1, B2, C1, C2 = (double3)(0.0);
    double F_left, F_right = 0.0;
    while(COUNTER < TERNARY_SEARCH_MAX_ITERATION && fabs(max - min) >= precision)
    {
        COUNTER += 1;
        min_third = min + (max - min) / 3.0;
        max_third = max - (max - min) / 3.0;

        shrink_tris(A, B, C, &A1, &B1, &C1, min_third, centroid_src, centroid_stretch);
        shrink_tris(A, B, C, &A2, &B2, &C2, max_third, centroid_src, centroid_stretch);
        F_left  = squared_error_tris(A1, B1, C1, QEM, idx);
        F_right = squared_error_tris(A2, B2, C2, QEM, idx);

        if(F_left>F_right) { min = min_third; }
        else               { max = max_third; }
    }
    
    return (min + max) / 2.0;
}

static inline double ternary_search_quad(double min, double max, double precision,
    const double3 A, const double3 B, const double3 C, const double3 D, // stretch quad
    double* QEM, const int idx, const double3 centroid_src, const double3 centroid_stretch)
{
    int COUNTER = 0;
    double min_third, max_third = 0.0;
    double3 A1, A2, B1, B2, C1, C2, D1, D2 = (double3)(0.0);
    double F_left, F_right = 0.0;
    while(COUNTER < TERNARY_SEARCH_MAX_ITERATION && fabs(max - min) >= precision)
    {
        COUNTER += 1;
        min_third = min + (max - min) / 3.0;
        max_third = max - (max - min) / 3.0;

        shrink_quad(A, B, C, D, &A1, &B1, &C1, &D1, min_third, centroid_src, centroid_stretch);
        shrink_quad(A, B, C, D, &A2, &B2, &C2, &D2, max_third, centroid_src, centroid_stretch);
        F_left  = squared_error_quad(A1, B1, C1, D1, QEM, idx);
        F_right = squared_error_quad(A2, B2, C2, D2, QEM, idx);

        if(F_left>F_right) { min = min_third; }
        else               { max = max_third; }
    }
    
    return (min + max) / 2.0;
}

static inline double ternary_search_ngon(double min, double max, double precision,
    double* P_stretch, double* P, double* P1, double* QEM,
    const int p_idx, const int q_idx, const int size,
    const double3 centroid_src, const double3 centroid_stretch)
{
    int COUNTER = 0;
    double min_third, max_third = 0.0;
    double F_left, F_right = 0.0;
    while(COUNTER < TERNARY_SEARCH_MAX_ITERATION && fabs(max - min) >= precision)
    {
        COUNTER += 1;
        min_third = min + (max - min) / 3.0;
        max_third = max - (max - min) / 3.0;

        shrink_ngon(P_stretch, P,  size, p_idx, min_third, centroid_src, centroid_stretch);
        shrink_ngon(P_stretch, P1, size, p_idx, max_third, centroid_src, centroid_stretch);
        F_left  = squared_error_ngon(P , QEM, p_idx, q_idx, size);
        F_right = squared_error_ngon(P1, QEM, p_idx, q_idx, size);

        if(F_left>F_right) { min = min_third; }
        else               { max = max_third; }
    }
    
    return (min + max) / 2.0;
}

#endif
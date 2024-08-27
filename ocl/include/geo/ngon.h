#ifndef __geo_ngon__
#define __geo_ngon__

/*  Library to extract geometric information from n-gon
    All methods are overloaded for float and double precision for convenience.

Methods:
    - ngon_quality: given an array of points, return the [0,1] quality value measuring how regular the ngon is. Optionally returns the perimeter as well
*/

#include "geo/tris.h"

static inline __attribute__((overloadable)) float ngon_quality(float* P_array, const int size, const int idx)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    float3 centroid = (float3)(0.0f);
    float q_max = 0.0f;
    float q_min = 1.0f;
    float q = 0.0f;
    for(int i=idx; i<idx+size; i++)
    {
        centroid += vload3(i, P_array);
    }
    centroid /= (float)(size);

    float d_min = 1.0f;
    float3 N_prev, N_curr = (float3)(0.0f);
    float3 A, B = (float3)(0.0f);
    for(int i=idx; i<idx+size; i++)
    {
        int i_next = i+1; if(i_next>=idx+size) { i_next = i_next - size; }
        A = vload3(i, P_array);
        B = vload3(i_next, P_array);
        q = tris_quality(A, B, centroid);
        if(q > q_max) { q_max = q; }
        if(q < q_min) { q_min = q; }

        if(i==0) // first iteration, just save current normal
        {
            N_prev = normalize(cross(A-centroid, B-centroid));
        }
        else // save maximum orientation delta between neighbouring triangles
        {
            N_curr = normalize(cross(A - centroid, B-centroid));
            d_min = fmin(d_min, clamp((float)(dot(N_prev, N_curr)), 0.0f, 1.0f));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

static inline __attribute__((overloadable)) double ngon_quality(double* P_array, const int size, const int idx)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    double3 centroid = (double3)(0.0);
    double q_max = 0.0;
    double q_min = 1.0;
    double q = 0.0;
    for(int i=idx; i<idx+size; i++)
    {
        centroid += vload3(i, P_array);
    }
    centroid /= (float)(size);

    double d_min = 1.0;
    double3 N_prev, N_curr = (double3)(0.0);
    double3 A, B = (double3)(0.0);
    for(int i=idx; i<idx+size; i++)
    {
        int i_next = i+1; if(i_next>=idx+size) { i_next = i_next - size; }
        A = vload3(i, P_array);
        B = vload3(i_next, P_array);
        q = tris_quality(A, B, centroid);
        if(q > q_max) { q_max = q; }
        if(q < q_min) { q_min = q; }

        if(i==0) // first iteration, just save current normal
        {
            N_prev = normalize(cross(A-centroid, B-centroid));
        }
        else // save maximum orientation delta between neighbouring triangles
        {
            N_curr = normalize(cross(A-centroid, B-centroid));
            d_min = fmin(d_min, clamp((double)(dot(N_prev, N_curr)), 0.0, 1.0));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

static inline __attribute__((overloadable)) float ngon_quality(float* P_array, const int size, const int idx, float* perimeter)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    float3 centroid = (float3)(0.0f);
    float q_max = 0.0f;
    float q_min = 1.0f;
    float q = 0.0f;
    for(int i=idx; i<idx+size; i++)
    {
        centroid += vload3(i, P_array);
    }
    centroid /= (float)(size);

    float d_min = 1.0;
    float3 N_prev, N_curr = (float3)(0.0f);
    float3 A, B = (float3)(0.0f);
    for(int i=idx; i<idx+size; i++)
    {
        int i_next = i+1; if(i_next>=idx+size) { i_next = i_next - size; }
        A = vload3(i, P_array);
        B = vload3(i_next, P_array);
        q = tris_quality(A, B, centroid);
        if(q > q_max) { q_max = q; }
        if(q < q_min) { q_min = q; }
        *perimeter += distance(A, B);

        if(i==0) // first iteration, just save current normal
        {
            N_prev = normalize(cross(A-centroid, B-centroid));
        }
        else // save maximum orientation delta between neighbouring triangles
        {
            N_curr = normalize(cross(A-centroid, B-centroid));
            d_min = fmin(d_min, clamp((float)(dot(N_prev, N_curr)), 0.0f, 1.0f));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

static inline __attribute__((overloadable)) double ngon_quality(double* P_array, const int size, const int idx, double* perimeter)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    double3 centroid = (double3)(0.0);
    double q_max = 0.0;
    double q_min = 1.0;
    double q = 0.0;
    for(int i=idx; i<idx+size; i++)
    {
        centroid += vload3(i, P_array);
    }
    centroid /= (float)(size);

    double d_min = 1.0;
    double3 N_prev, N_curr = (double3)(0.0);
    double3 A, B = (double3)(0.0);
    for(int i=idx; i<idx+size; i++)
    {
        int i_next = i+1; if(i_next>=idx+size) { i_next = i_next - size; }
        A = vload3(i, P_array);
        B = vload3(i_next, P_array);
        q = tris_quality(A, B, centroid);
        if(q > q_max) { q_max = q; }
        if(q < q_min) { q_min = q; }
        *perimeter += distance(A, B);

        if(i==0) // first iteration, just save current normal
        {
            N_prev = normalize(cross(A-centroid, B-centroid));
        }
        else // save maximum orientation delta between neighbouring triangles
        {
            N_curr = normalize(cross(A-centroid, B-centroid));
            d_min = fmin(d_min, clamp((double)(dot(N_prev, N_curr)), 0.0, 1.0));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

#endif
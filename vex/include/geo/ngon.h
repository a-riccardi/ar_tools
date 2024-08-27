#ifndef __geo_ngon__
#define __geo_ngon__

/*  Library to extract geometric information from n-gon

Methods:
    - ngon_quality: given an array of points, return the [0,1] quality value measuring how regular the ngon is. Optionally returns the perimeter and/or centroid
*/

#include "geo/tris.h"

function float ngon_quality(vector P_array[]; int size)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    vector centroid = (vector)(0.0);
    float q_max = 0.0;
    float q_min = 1.0;
    float q = 0.0;
    for(int i=0; i<size; i++)
    {
        centroid += P_array[i];
    }
    centroid /= (float)(size);

    float d_min = 1.0;
    vector N_prev, N_curr = (vector)(0.0);
    vector A, B = (vector)(0.0);
    for(int i=0; i<size; i++)
    {
        int i_next = i+1; if(i_next>=size) { i_next = i_next - size; }
        A = P_array[i];
        B = P_array[i_next];
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
            d_min = min(d_min, clamp(dot(N_prev, N_curr), 0.0, 1.0));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

function float ngon_quality(vector P_array[]; int size; export vector centroid)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    float q_max = 0.0;
    float q_min = 1.0;
    float q = 0.0;
    for(int i=0; i<size; i++)
    {
        centroid += P_array[i];
    }
    centroid /= (float)(size);

    float d_min = 1.0;
    vector N_prev, N_curr = (vector)(0.0);
    vector A, B = (vector)(0.0);
    for(int i=0; i<size; i++)
    {
        int i_next = i+1; if(i_next>=size) { i_next = i_next - size; }
        A = P_array[i];
        B = P_array[i_next];
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
            d_min = min(d_min, clamp(dot(N_prev, N_curr), 0.0, 1.0));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

function float ngon_quality(vector P_array[]; int size; export float perimeter)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    vector centroid = (vector)(0.0);
    float q_max = 0.0;
    float q_min = 1.0;
    float q = 0.0;
    for(int i=0; i<size; i++)
    {
        centroid += P_array[i];
    }
    centroid /= (float)(size);

    float d_min = 1.0;
    vector N_prev, N_curr = (vector)(0.0);
    vector A, B = (vector)(0.0);
    for(int i=0; i<size; i++)
    {
        int i_next = i+1; if(i_next>=size) { i_next = i_next - size; }
        A = P_array[i];
        B = P_array[i_next];
        q = tris_quality(A, B, centroid);
        if(q > q_max) { q_max = q; }
        if(q < q_min) { q_min = q; }
        perimeter += distance(A, B);

        if(i==0) // first iteration, just save current normal
        {
            N_prev = normalize(cross(A-centroid, B-centroid));
        }
        else // save maximum orientation delta between neighbouring triangles
        {
            N_curr = normalize(cross(A-centroid, B-centroid));
            d_min = min(d_min, clamp(dot(N_prev, N_curr), 0.0, 1.0));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

function float ngon_quality(vector P_array[]; int size; export float perimeter; export vector centroid)
{
    /* Measure the quality of the n-gon. After building the ngon centroid we measure every triangle P-Pnext-Centroid quality,
    returning the ratio between the worst and best quality. It also checks how much coplanar the triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how regular the ngon is.
    */

    float q_max = 0.0;
    float q_min = 1.0;
    float q = 0.0;
    for(int i=0; i<size; i++)
    {
        centroid += P_array[i];
    }
    centroid /= (float)(size);

    float d_min = 1.0;
    vector N_prev, N_curr = (vector)(0.0);
    vector A, B = (vector)(0.0);
    for(int i=0; i<size; i++)
    {
        int i_next = i+1; if(i_next>=size) { i_next = i_next - size; }
        A = P_array[i];
        B = P_array[i_next];
        q = tris_quality(A, B, centroid);
        if(q > q_max) { q_max = q; }
        if(q < q_min) { q_min = q; }
        perimeter += distance(A, B);

        if(i==0) // first iteration, just save current normal
        {
            N_prev = normalize(cross(A-centroid, B-centroid));
        }
        else // save maximum orientation delta between neighbouring triangles
        {
            N_curr = normalize(cross(A-centroid, B-centroid));
            d_min = min(d_min, clamp(dot(N_prev, N_curr), 0.0, 1.0));
            N_prev = N_curr;
        }
    }

    return (q_min / q_max) * d_min;
}

#endif
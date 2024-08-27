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
                           Overload available for a QEM in matrix form.
    - squared_error_tris: returns the total squared distance error for each point of the triangle.
    - squared_error_quad: returns the total squared distance error for each point of the quad.
    - squared_error_ngon: returns the total squared distance error for each point of the n-gon.
    - ternary_search_tris: returns the optimal value of `k` which minimizes the squared error distance
    - ternary_search_quad: returns the optimal value of `k` which minimizes the squared error distance
    - ternary_search_ngon: returns the optimal value of `k` which minimizes the squared error distance
*/

#include "geo/tracing.h"

#define TERNARY_SEARCH_MAX_ITERATION 100

function void project_on_plane_tris(export vector A; export vector B; export vector C;
   vector prim_centroid; vector prim_normal)
{
    A = project_on_plane(A, prim_centroid, prim_normal);
    B = project_on_plane(B, prim_centroid, prim_normal);
    C = project_on_plane(C, prim_centroid, prim_normal);
}

function void project_on_plane_tris(export vector A; export vector B; export vector C;
   vector4 abcd)
{
    A = project_on_plane(A, abcd);
    B = project_on_plane(B, abcd);
    C = project_on_plane(C, abcd);
}

function void project_on_plane_quad(export vector A; export vector B; export vector C; export vector D;
   vector prim_centroid; vector prim_normal)
{
    A = project_on_plane(A, prim_centroid, prim_normal); 
    B = project_on_plane(B, prim_centroid, prim_normal);
    C = project_on_plane(C, prim_centroid, prim_normal);
    D = project_on_plane(D, prim_centroid, prim_normal);
}

function void project_on_plane_quad(export vector A; export vector B; export vector C; export vector D;
   vector4 abcd)
{
    A = project_on_plane(A, abcd); 
    B = project_on_plane(B, abcd);
    C = project_on_plane(C, abcd);
    D = project_on_plane(D, abcd);
}

function void project_on_plane_ngon(vector P[]; vector P_planar[];
   vector prim_centroid; vector prim_normal)
{
    foreach(int i; vector P_src; P)
    {
        P_planar[i] = project_on_plane(P_src, prim_centroid, prim_normal);
    }
}

function void project_on_plane_ngon(vector P[]; vector P_planar[];
   vector4 abcd)
{
    foreach(int i; vector P_src; P)
    {
        P_planar[i] = project_on_plane(P_src, abcd);
    }
}

function void stretch_tris(vector A; vector B; vector C;
    export vector A_stretch; export vector B_stretch; export vector C_stretch;
    export vector centroid_src; export vector centroid_stretch; vector prim_normal; float lambda)
{
    centroid_src = (A + B + C) / 3.0;

    vector N  = (vector)(0.0);
    vector edge = (vector)(0.0);

    // extrude A
    edge = normalize(B-C);
    N = (cross(prim_normal, edge));
    A_stretch = A + N * lambda;

    // extrude B
    edge = normalize(C-A);
    N = (cross(prim_normal, edge));
    B_stretch = B + N * lambda;

    // extrude C
    edge = normalize(A-B);
    N = (cross(prim_normal, edge));
    C_stretch = C + N * lambda;

    centroid_stretch = (A_stretch + B_stretch + C_stretch) / 3.0;
}

function void stretch_quad(vector A; vector B; vector C; vector D;
    export vector A_stretch; export vector B_stretch; export vector C_stretch; export vector D_stretch;
    export vector centroid_src; export vector centroid_stretch; vector prim_normal; float lambda)
{
    centroid_src = (A + B + C + D) / 4.0;

    vector N  = (vector)(0.0);
    vector edge = (vector)(0.0);

    // extrude A
    edge = normalize(B-D);
    N = (cross(prim_normal, edge));
    A_stretch = A + N * lambda;

    // extrude B
    edge = normalize(C-A);
    N = (cross(prim_normal, edge));
    B_stretch = B + N * lambda;

    // extrude C
    edge = normalize(D-B);
    N = (cross(prim_normal, edge));
    C_stretch = C + N * lambda;

    // extrude D
    edge = normalize(A-C);
    N = (cross(prim_normal, edge));
    D_stretch = D + N * lambda;

    centroid_stretch = (A_stretch + B_stretch + C_stretch + D_stretch) / 4.0;
}

function void stretch_ngon(vector P[]; vector P_stretch[]; int size;
    export vector centroid_src; export vector centroid_stretch; vector prim_normal; float lambda)
{
    vector N = (vector)(0.0);
    vector edge = (vector)(0.0);
    vector P1 = (vector)(0.0);

    // re-initializing due to previous data inside? 
    centroid_src = (vector)(0.0);
    centroid_stretch = (vector)(0.0);

    int pt_curr, pt_prev, pt_next = 0;
    for(int i=0; i<size; i++)
    {
        // grab previous/current/next point indexes.
        pt_curr = i;
        pt_prev = i-1; if(pt_prev<0)     { pt_prev += size; }
        pt_next = i+1; if(pt_next>=size) { pt_next -= size; }

        edge = normalize(P[pt_next] - P[pt_prev]);
        N = cross(prim_normal, edge);
        P1 = P[pt_curr] + N * lambda;

        centroid_src += P[pt_curr];
        centroid_stretch += P1;
        P_stretch[pt_curr] = P1;
    }

    centroid_src /= (float)(size);
    centroid_stretch /= (float)(size);
}

function void shrink_tris(vector A; vector B; vector C;
    export vector A_shrink; export vector B_shrink; export vector C_shrink;
    float k; vector centroid_src; vector centroid_stretch)
{
    A_shrink = centroid_src + (A - centroid_stretch) * k;
    B_shrink = centroid_src + (B - centroid_stretch) * k;
    C_shrink = centroid_src + (C - centroid_stretch) * k;
}

function void shrink_quad(vector A; vector B; vector C; vector D;
    export vector A_shrink; export vector B_shrink; export vector C_shrink; export vector D_shrink;
    float k; vector centroid_src; vector centroid_stretch)
{
    A_shrink = centroid_src + (A - centroid_stretch) * k;
    B_shrink = centroid_src + (B - centroid_stretch) * k;
    C_shrink = centroid_src + (C - centroid_stretch) * k;
    D_shrink = centroid_src + (D - centroid_stretch) * k;
}

function void shrink_ngon(vector P_stretch[]; vector P_shrink[];
    float k; vector centroid_src; vector centroid_stretch)
{
    foreach(int i; vector P; P_stretch)
    {
        P_shrink[i] = centroid_src + (P - centroid_stretch) * k;
    }
}

function float pt_surface_distance(vector P; float QEM[]; int offset)
{
    /* Given a position P and the quadric error matrix QEM, computes the quadric error */

    return QEM[offset+0]*P.x*P.x + 2.0*QEM[offset+1]*P.x*P.y + 2.0*QEM[offset+2]*P.x*P.z + 2.0*QEM[offset+3]*P.x
                                 +     QEM[offset+4]*P.y*P.y + 2.0*QEM[offset+5]*P.y*P.z + 2.0*QEM[offset+6]*P.y
                                                             +     QEM[offset+7]*P.z*P.z + 2.0*QEM[offset+8]*P.z 
                                                                                           +   QEM[offset+9];
}

function float pt_surface_distance(vector P; matrix QEM)
{
    /* Given a position P and the quadric error matrix QEM, computes the quadric error */

    return QEM.xx*P.x*P.x + 2.0*QEM.xy*P.x*P.y + 2.0*QEM.xz*P.x*P.z + 2.0*QEM.xw*P.x
                          +     QEM.yy*P.y*P.y + 2.0*QEM.yz*P.y*P.z + 2.0*QEM.yw*P.y
                                               +     QEM.zz*P.z*P.z + 2.0*QEM.zw*P.z 
                                                                    +     QEM.ww;
}

function float squared_error_tris(vector A; vector B; vector C; float QEM[])
{
    return pt_surface_distance(A, QEM, 0 )
         + pt_surface_distance(B, QEM, 10)
         + pt_surface_distance(C, QEM, 20);

}

function float squared_error_quad(vector A; vector B; vector C; vector D; float QEM[])
{
    return pt_surface_distance(A, QEM, 0 )
         + pt_surface_distance(B, QEM, 10)
         + pt_surface_distance(C, QEM, 20)
         + pt_surface_distance(D, QEM, 30);
}

function float squared_error_ngon(vector P[]; float QEM[])
{
    float error = 0.0;
    foreach(int i; vector P_offset; P)
    {
        error += pt_surface_distance(P_offset, QEM, i*10);
    }
    return error;
}

function float ternary_search_tris(float min; float max; float precision;
    vector A; vector B; vector C; // stretch triangle
    float QEM[]; vector centroid_src; vector centroid_stretch)
{
    int COUNTER = 0;
    float min_third, max_third = 0.0;
    vector A1, A2, B1, B2, C1, C2 = (vector)(0.0);
    float F_left, F_right = 0.0;
    while(COUNTER < TERNARY_SEARCH_MAX_ITERATION && abs(max - min) >= precision)
    {
        COUNTER += 1;
        min_third = min + (max - min) / 3.0;
        max_third = max - (max - min) / 3.0;

        shrink_tris(A, B, C, A1, B1, C1, min_third, centroid_src, centroid_stretch);
        shrink_tris(A, B, C, A2, B2, C2, max_third, centroid_src, centroid_stretch);
        F_left  = squared_error_tris(A1, B1, C1, QEM);
        F_right = squared_error_tris(A2, B2, C2, QEM);

        if(F_left>F_right) { min = min_third; }
        else               { max = max_third; }
    }
    
    return (min + max) / 2.0;
}

function float ternary_search_quad(float min; float max; float precision;
    vector A; vector B; vector C; vector D; // stretch quad
    float QEM[]; vector centroid_src; vector centroid_stretch)
{
    int COUNTER = 0;
    float min_third, max_third = 0.0;
    vector A1, A2, B1, B2, C1, C2, D1, D2 = (vector)(0.0);
    float F_left, F_right = 0.0;
    while(COUNTER < TERNARY_SEARCH_MAX_ITERATION && abs(max - min) >= precision)
    {
        COUNTER += 1;
        min_third = min + (max - min) / 3.0;
        max_third = max - (max - min) / 3.0;

        shrink_quad(A, B, C, D, A1, B1, C1, D1, min_third, centroid_src, centroid_stretch);
        shrink_quad(A, B, C, D, A2, B2, C2, D2, max_third, centroid_src, centroid_stretch);
        F_left  = squared_error_quad(A1, B1, C1, D1, QEM);
        F_right = squared_error_quad(A2, B2, C2, D2, QEM);

        if(F_left>F_right) { min = min_third; }
        else               { max = max_third; }
    }
    
    return (min + max) / 2.0;
}

function float ternary_search_ngon(float min; float max; float precision;
    vector P_stretch[]; vector P[]; vector P1[]; float QEM[];
    vector centroid_src; vector centroid_stretch)
{
    int COUNTER = 0;
    float min_third, max_third = 0.0;
    float F_left, F_right = 0.0;
    while(COUNTER < TERNARY_SEARCH_MAX_ITERATION && abs(max - min) >= precision)
    {
        COUNTER += 1;
        min_third = min + (max - min) / 3.0;
        max_third = max - (max - min) / 3.0;

        shrink_ngon(P_stretch, P,  min_third, centroid_src, centroid_stretch);
        shrink_ngon(P_stretch, P1, max_third, centroid_src, centroid_stretch);
        F_left  = squared_error_ngon(P , QEM);
        F_right = squared_error_ngon(P1, QEM);

        if(F_left>F_right) { min = min_third; }
        else               { max = max_third; }
    }
    
    return (min + max) / 2.0;
}

#endif
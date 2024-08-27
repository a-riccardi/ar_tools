#include "geo/tris.h"
#include "geo/quad.h"
#include "geo/ngon.h"
#include "geo/gem_sso.h"

#define K_SEARCH_PRECISION 0.0000001
#define TERNARY_SEARCH_MAX_ITERATION 25

kernel void GED_SSO( 
#ifdef SINGLE_WORKGROUP
#ifdef SINGLE_WORKGROUP_SPANS
    int cur_workset,
#endif
    int num_worksets,
    global const int *worksets_begins,
    global const int *worksets_lengths,
#else
    int worksets_begin,
    int worksets_length,
#endif
    int mode,
    float treshold,
    int primpoints_length,
    global int* restrict primpoints_index,
    global int* restrict primpoints,
    int quality_length,
    global double* restrict quality,
    int quality_delta_length,
    global double* restrict quality_delta,
    int shrink_factor_length,
    global double* restrict shrink_factor,
    int Plane_length,
    global double* restrict Plane,
#ifdef HAS_P_planar
    int P_planar_length,
    global int* restrict P_planar_index,
    global double* restrict P_planar,
#endif
#ifdef HAS_P_stretch
    int P_stretch_length,
    global int* restrict P_stretch_index,
    global double* restrict P_stretch,
#endif
#ifdef HAS_P_buffer
    int P_buffer_length,
    global int* restrict P_buffer_index,
    global double* restrict P_buffer,
#endif
    int P_shrink_length,
    global int* restrict P_shrink_index,
    global double* restrict P_shrink,
    int QEM_length,
    global int* restrict QEM_index,
    global double* restrict QEM,
    int Pos_length,
    global double* restrict Pos,
    int is_unshared_length,
    global int * restrict is_unshared,
    int pointprims_length,
    global int* restrict pointprims_index,
    global int* restrict pointprims,
    int primpoint_idx_length,
    global int* restrict primpoint_idx_index,
    global int* restrict primpoint_idx
)
{
    // iterate on each primitive to compute the new Stretch-Shrink Operation position
    // for each of the vertices

#ifdef SINGLE_WORKGROUP
#ifdef SINGLE_WORKGROUP_SPANS
   for(int i=cur_workset; i<cur_workset+num_worksets; i++)
#else
   for(int i=0; i<num_worksets; i++)
#endif
    {
        int worksets_begin = worksets_begins[i];
        int worksets_length = worksets_lengths[i];
        if (i > 0) { barrier(CLK_GLOBAL_MEM_FENCE); }
#else
    { 
#endif
        int idx = get_global_id(0);
        if (idx < worksets_length) { idx += worksets_begin; }
        if (idx >= quality_length) { return;                }

        int q_start = QEM_index[idx];
        int n_start = primpoints_index[idx];
        int n_end = primpoints_index[idx+1];
        int size = n_end - n_start;
        double fsize = (double)(size);

        // l: perimeter
        // q: triangle quality
        // q_delta = 1.0 + delta q
        double l = 0.0;
        double q = 1.0;

        double4 abcd = vload4(idx, Plane);

        if(size == 3)
        {
            // current primitive is a triangle
            double3 A = vload3(primpoints[n_start+0], Pos);
            double3 B = vload3(primpoints[n_start+1], Pos);
            double3 C = vload3(primpoints[n_start+2], Pos);
            
            project_on_plane_tris(&A, &B, &C, abcd);
            
            // get current quality and perimeter
            q = tris_quality(A, B, C, &l);
            l /= 3.0;

            if(q>=treshold)
            {
                // if the triangle is already sufficiently equilateral
                // just store the current vertex positions in P_shrink 
                vstore3(A, n_start+0, P_shrink);
                vstore3(B, n_start+1, P_shrink);
                vstore3(C, n_start+2, P_shrink);

                quality[idx] = q;
                quality_delta[idx] = 0.0;
                shrink_factor[idx] = 1.0;
            }
            else
            {
                // compute stretch factor
                double lambda = (1.0 - q) * l;
                double3 A1, B1, C1 = (double3)(0.0);
                double3 centroid_src, centroid_stretch = (double3)(0.0);
                stretch_tris(A, B, C, &A1, &B1, &C1, &centroid_src, &centroid_stretch, abcd.xyz, lambda);

                double l_stretch, q_stretch = 0.0;
                q_stretch = tris_quality(A1, B1, C1, &l_stretch);
                l_stretch /= 3.0;

                // k is the shrinking ratio after stretching
                // pre-initializing k to the default planar patch case
                double k = l / l_stretch;

                double3 A_shrink, B_shrink, C_shrink = (double3)(0.0);

                if(mode == 1)
                {
                    // mesh surface case:
                    double min =-2.0;
                    double max = 2.0;

                    k = ternary_search_tris(min, max, K_SEARCH_PRECISION, A1, B1, C1, QEM, q_start, centroid_src, centroid_stretch);
                }

                shrink_tris(A1, B1, C1, &A_shrink, &B_shrink, &C_shrink, k, centroid_src, centroid_stretch);

                vstore3(A_shrink, n_start+0, P_shrink);
                vstore3(B_shrink, n_start+1, P_shrink);
                vstore3(C_shrink, n_start+2, P_shrink);

                shrink_factor[idx] = k;
                quality_delta[idx] = q_stretch - quality[idx];
                quality[idx] = q;
            }
        }
        else if(size == 4)
        {
            // current primitive is a quad
            double3 A = vload3(primpoints[n_start+0], Pos);
            double3 B = vload3(primpoints[n_start+1], Pos);
            double3 C = vload3(primpoints[n_start+2], Pos);
            double3 D = vload3(primpoints[n_start+3], Pos);

            project_on_plane_quad(&A, &B, &C, &D, abcd);

            // get current quality and perimeter
            q = quad_quality(A, B, C, D, &l);
            l /= 4.0;

            if(q>=treshold)
            {
                // if the quad is already sufficiently equilateral
                // just store the current vertex positions in P_shrink
                vstore3(A, n_start+0, P_shrink);
                vstore3(B, n_start+1, P_shrink);
                vstore3(C, n_start+2, P_shrink);
                vstore3(D, n_start+3, P_shrink);

                quality[idx] = q;
                quality_delta[idx] = 0.0;
                shrink_factor[idx] = 1.0;
            }
            else
            {
                // compute stretch factor
                double lambda = (1.0 - q) * l;
                double3 A1, B1, C1, D1 = (double3)(0.0);
                double3 centroid_src, centroid_stretch = (double3)(0.0);
                stretch_quad(A, B, C, D, &A1, &B1, &C1, &D1, &centroid_src, &centroid_stretch, abcd.xyz, lambda);

                double l_stretch, q_stretch = 0.0;
                q_stretch = quad_quality(A1, B1, C1, D1, &l_stretch);
                l_stretch /= 4.0;
                // k is the shrinking ratio after stretching
                // pre-initializing k to the default planar patch case
                double k = l / l_stretch;

                double3 A_shrink, B_shrink, C_shrink, D_shrink = (double3)(0.0);
                
                if(mode == 1)
                {
                    // mesh surface case:
                    double min =-2.0;
                    double max = 2.0;

                    k = ternary_search_quad(min, max, K_SEARCH_PRECISION, A1, B1, C1, D1, QEM, q_start, centroid_src, centroid_stretch);
                }

                shrink_quad(A1, B1, C1, D1, &A_shrink, &B_shrink, &C_shrink, &D_shrink, k, centroid_src, centroid_stretch);

                vstore3(A_shrink, n_start+0, P_shrink);
                vstore3(B_shrink, n_start+1, P_shrink);
                vstore3(C_shrink, n_start+2, P_shrink);
                vstore3(D_shrink, n_start+3, P_shrink);
                
                shrink_factor[idx] = k;
                quality_delta[idx] = q_stretch - quality[idx];
                quality[idx] = q;
            }
        }
        else
        {
#ifdef HAS_P_planar
            // current polygon is an n-gon
            project_on_plane_ngon(primpoints, Pos, P_planar, size, n_start, abcd);
            
            // get current quality and perimeter
            q = ngon_quality(P_planar, size, n_start, &l);
            l /= fsize;

            if(q>=treshold)
            {
                // if the ngon is already sufficiently equilateral
                // just store the current vertex positions in P_shrink
                for(int i=n_start; i<n_start+size; i++)
                {
                    vstore3(vload3(i, P_planar), i, P_shrink);
                }

                quality[idx] = q;
                quality_delta[idx] = 0.0;
                shrink_factor[idx] = 1.0;
            }
            else
            {
                // compute stretch factor
                double lambda = (1.0 - q) * l;
                
                double3 centroid_src, centroid_stretch = (double3)(0.0);
                stretch_ngon(P_planar, P_stretch, size, n_start, &centroid_src, &centroid_stretch, abcd.xyz, lambda);

                double l_stretch, q_stretch = 0.0;
                q_stretch = ngon_quality(P_stretch, size, n_start, &l_stretch);
                l_stretch /= fsize;

                // k is the shrinking ratio after stretching
                // pre-initializing k to the default planar patch case
                double k = l / l_stretch;

                if(mode == 1)
                {
                    // mesh surface case:
                    double min =-2.0;
                    double max = 2.0;

                    k = ternary_search_ngon(min, max, K_SEARCH_PRECISION, P_stretch, P_shrink, P_buffer, QEM, n_start, q_start, size, centroid_src, centroid_stretch);
                }

                shrink_ngon(P_stretch, P_shrink, size, n_start, k, centroid_src, centroid_stretch);

                shrink_factor[idx] = k;
                quality_delta[idx] = q_stretch - quality[idx];
                quality[idx] = q;
            }
#else
            for(int i=n_start; i<n_start+size; i++)
            {
                vstore3(vload3(i, Pos), i, P_shrink);
            }

            quality[idx] = q;
            quality_delta[idx] = 0.0;
            shrink_factor[idx] = 1.0;
#endif
        }
    }    
}
 
kernel void GED_position_weight( 
#ifdef SINGLE_WORKGROUP
#ifdef SINGLE_WORKGROUP_SPANS
    int cur_workset,
#endif
    int num_worksets,
    global const int *worksets_begins,
    global const int *worksets_lengths,
#else
    int worksets_begin,
    int worksets_length,
#endif
    int mode,
    float treshold,
    int primpoints_length,
    global int* restrict primpoints_index,
    global int* restrict primpoints,
    int quality_length,
    global double* restrict quality,
    int quality_delta_length,
    global double* restrict quality_delta,
    int shrink_factor_length,
    global double* restrict shrink_factor,
    int Plane_length,
    global double* restrict Plane,
#ifdef HAS_P_planar
    int P_planar_length,
    global int* restrict P_planar_index,
    global double* restrict P_planar,
#endif
#ifdef HAS_P_stretch
    int P_stretch_length,
    global int* restrict P_stretch_index,
    global double* restrict P_stretch,
#endif
#ifdef HAS_P_buffer
    int P_buffer_length,
    global int* restrict P_buffer_index,
    global double* restrict P_buffer,
#endif
    int P_shrink_length,
    global int* restrict P_shrink_index,
    global double* restrict P_shrink,
    int QEM_length,
    global int* restrict QEM_index,
    global double* restrict QEM,
    int Pos_length,
    global double* restrict Pos,
    int is_unshared_length,
    global int * restrict is_unshared,
    int pointprims_length,
    global int* restrict pointprims_index,
    global int* restrict pointprims,
    int primpoint_idx_length,
    global int* restrict primpoint_idx_index,
    global int* restrict primpoint_idx
)
{
    // iterate on each point to compute the weighted average of the 
    // primitives new position
#ifdef SINGLE_WORKGROUP
#ifdef SINGLE_WORKGROUP_SPANS
   for(int i=cur_workset; i<cur_workset+num_worksets; i++)
#else
   for(int i=0; i<num_worksets; i++)
#endif
    {
        int worksets_begin = worksets_begins[i];
        int worksets_length = worksets_lengths[i];
        if (i > 0) { barrier(CLK_GLOBAL_MEM_FENCE); }
#else
    {
#endif
        int idx = get_global_id(0);
        if (idx < worksets_length) { idx += worksets_begin; }
        if (idx >= Pos_length)     { return;                } 
        if (is_unshared[idx])      { return;                } 

        double3 P = (double3)(0.0);
        double weight = 0.0;

        int prim_start = pointprims_index[idx];
        int prim_end = pointprims_index[idx+1];
        for(int i=prim_start; i<prim_end; i++)
        {
            int prim = pointprims[i];
            int point_idx = primpoint_idx[i];
            int prim_offset = primpoints_index[prim]; 

            double q_delta = 1.0 + quality_delta[prim];
            
            P += vload3(prim_offset + point_idx, P_shrink) * q_delta;
            weight += q_delta;    
        }
        
        vstore3(P / weight, idx, Pos);
    }
}

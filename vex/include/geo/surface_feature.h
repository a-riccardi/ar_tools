#ifndef __geo_surface_feature__
#define __geo_surface_feature__

/* This module contains functions to analize surface features like occlusion and thickness

Methods:
    - trace_ws: send rays from `P_src` along all provided `directions[]` - it expects the direction array to be in world space
    - trace_ts: send rays from `P_src` along all provided `directions[]` - it expects the direction array to be in tangent space
    - trace_occlusion_thickness: send rays from `P_src` both outward and inward along all provided `directions[]`.
                                 If `rtype` is 0 returns the occlusion/thickness [0,1] value, otherwise returns the minimum hit distance
*/


function float trace_ws(
    const vector P_src;
    const vector offset;
    const vector directions[];
    const float max_distance;
    const int rtype
)
{
    /* Traces a set of rays defined by the directions array, testing self-intersection.
    Depending on the rtype value, it either provides a coverage % value in [0..1] range
    or the minimum hit distance.
    
    - P_src: the ray start position in world space
    - directions: a world-space array of directions to cast rays to
    - offset: a world-space offset to apply to P_src along distance to ensure no self-intersection
    - max_distance: tracing distance for the intersection test
    - rtype: what to return from the function, one of:
        - 0: return the coverage % in [0..1] range
        - 1: return the minimum distance hit
    */
    
    if(len(directions) < 1) { return max_distance; }

    float _search_distance = max_distance * (rtype == 0 ? 1.0f : 2.0f);
    float _min_distance = _search_distance * _search_distance;
    float _hit_counter = 0.0f;
    vector P_hit, uvw;
    foreach(vector ray_dir; directions)
    {
        vector ray = _search_distance * ray_dir;
        
        if(intersect(0, P_src + offset, ray, P_hit, uvw) > -1)
        {
            _min_distance = min(_min_distance, distance2(P_src, P_hit));
            _hit_counter += 1.0f;
        }
    }

    return rtype == 0 ? 
        (len(directions) - _hit_counter) / (float)(len(directions)) :
        clamp(sqrt(_min_distance), 0.0f, max_distance) * (_hit_counter > 0 ? 0.5f : 1.0f);
}

function float trace_ts(
    const vector P_src;
    const vector offset;
    const vector directions[];
    const vector4 orient;
    const float max_distance;
    const int rtype
)
{
    /* Traces a set of rays defined by the directions array, testing self-intersection.
    Depending on the rtype value, it either provides a coverage % value in [0..1] range
    or the minimum hit distance.

    - P_src: the ray start position in world space
    - directions: a tangent-space array of directions to cast rays to, Z-up
    - orient: a quaternion attribute to rotate the tangent-space directions to world-space 
    - offset: a world-space offset to apply to P_src along distance to ensure no self-intersection
    - max_distance: tracing distance for the intersection test
    - rtype: what to return from the function, one of:
        - 0: return the coverage % in [0..1] range
        - 1: return the minimum distance hit
    */

    if(len(directions) < 1) { return max_distance; }

    float _search_distance = max_distance * (rtype == 0 ? 1.0f : 2.0f);
    float _min_distance = _search_distance * _search_distance;
    float _hit_counter = 0.0f;
    vector P_hit, uvw;
    foreach(vector ray_dir; directions)
    {
        vector ray = _search_distance * qrotate(orient, ray_dir);
        
        if(intersect(0, P_src + offset, ray, P_hit, uvw) > -1)
        {
            _min_distance = min(_min_distance, distance2(P_src, P_hit));
            _hit_counter += 1.0f;
        }
    }

    return rtype == 0 ? 
        (len(directions) - _hit_counter) / (float)(len(directions)) :
        clamp(sqrt(_min_distance), 0.0f, abs(max_distance)) * (_hit_counter > 0 ? 0.5f : 1.0f);
}

function void trace_occlusion_thickness(
    const vector P_src;
    const vector offset;
    const vector directions[];
    const vector4 orient;
    const float max_distance;
    const int rtype;
    export float occlusion;
    export float thickness
)
{
    /* Traces a set of rays defined by the directions array both forward and backward,
    testing occlusion and thickness at the same time.
    Depending on the rtype value, it either provides a coverage % value in [0..1] range
    or the minimum hit distance.

    - P_src: the ray start position in world space
    - directions: a tangent-space array of directions to cast rays to, Z-up
    - orient: a quaternion attribute to rotate the tangent-space directions to world-space 
    - offset: a world-space offset to apply to P_src along distance to ensure no self-intersection
    - max_distance: tracing distance for the intersection test
    - rtype: what to return from the function, one of:
        - 0: return the coverage % in [0..1] range
        - 1: return the minimum distance hit
    */

    if(len(directions) < 1)
    {
        occlusion = max_distance;
        thickness = max_distance;
        return;
    }
    
    float _search_distance = max_distance * (rtype == 0 ? 1.0f : 2.0f);
    float _occlusion = _search_distance * _search_distance;
    float _thickness = _search_distance * _search_distance;
    float occlusion_hit = 0.0f;
    float thickness_hit = 0.0f;
    vector P_hit, uvw;

    for(int i=0; i<len(directions); i++)
    {
        vector ray = _search_distance * qrotate(orient, directions[i]);
        
        if(intersect(0, P_src + offset, ray, P_hit, uvw) > -1)
        {
            _occlusion = min(_occlusion, distance2(P_src, P_hit));
            occlusion_hit += 1.0f;
        }
            
        if(intersect(0, P_src - offset, -ray, P_hit, uvw) > -1)
        {
            _thickness = min(_thickness, distance2(P_src, P_hit));
            thickness_hit += 1.0f;
        }
    }
    
    occlusion = rtype == 0 ? 
        (occlusion_hit) / (float)(len(directions)) :
        clamp(sqrt(_occlusion), 0.0f, max_distance) * (occlusion_hit > 0 ? 0.5f : 1.0f);
    thickness = rtype == 0 ?
        (len(directions) - thickness_hit) / (float)(len(directions)) :
        clamp(sqrt(_thickness), 0.0f, max_distance) * (thickness_hit > 0 ? 0.5f : 1.0f);
}

#endif
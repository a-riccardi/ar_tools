#ifndef __geo_prims__
#define __geo_prims__

/*  Library to handle primitive attributes retrieval and iterpolation, mirroring vex functionality.
    All methods are overloaded for float and double precision for convenience.

Methods:
    - primuv_float:  given a point float attribute array, a primitive number and a uvw barycentric coordinates, interpolate the point attribute values
                     at the parametric coords. Works for tris, quads and n-gon
    - primuv_float2: given a point float2 attribute array, a primitive number and a uvw barycentric coordinates, interpolate the point attribute values
                     at the parametric coords. Works for tris, quads and n-gon
    - primuv_float3: given a point float3 attribute array, a primitive number and a uvw barycentric coordinates, interpolate the point attribute values
                     at the parametric coords. Works for tris, quads and n-gon
    - primuv_float4: given a point float4 attribute array, a primitive number and a uvw barycentric coordinates, interpolate the point attribute values
                     at the parametric coords. Works for tris, quads and n-gon
    - prim_float:  given a point float attribute array and a primitive number, interpolate the point attribute values in the center of the primitive - 
                   equivalent as primuv_float with {0.5, 0.5, 0.0} uvw coords. Works for tris, quads and n-gon
    - prim_float2: given a point float2 attribute array and a primitive number, interpolate the point attribute values in the center of the primitive - 
                   equivalent as primuv_float2 with {0.5, 0.5, 0.0} uvw coords. Works for tris, quads and n-gon
    - prim_float3: given a point float3 attribute array and a primitive number, interpolate the point attribute values in the center of the primitive - 
                   equivalent as primuv_float3 with {0.5, 0.5, 0.0} uvw coords. Works for tris, quads and n-gon
    - prim_float4: given a point float4 attribute array and a primitive number, interpolate the point attribute values in the center of the primitive - 
                   equivalent as primuv_float4 with {0.5, 0.5, 0.0} uvw coords. Works for tris, quads and n-gon
*/

static inline __attribute__((overloadable)) float primuv_float(
    const int geo, global float* attr, const int primnum, const float3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float _attr = 0.0f;
    float size = (float)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = attr[primpoints[pt_str  ]] * (1.0f - uvw.x - uvw.y) 
              + attr[primpoints[pt_str+1]] * uvw.x
              + attr[primpoints[pt_str+2]] * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = attr[primpoints[pt_str  ]] * (1.0f - uvw.x) * (1.0f - uvw.y) 
              + attr[primpoints[pt_str+1]] * (1.0f - uvw.x) * uvw.y
              + attr[primpoints[pt_str+2]] * uvw.x          * uvw.y
              + attr[primpoints[pt_str+3]] * uvw.x          * (1.0f - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        float v_accum_scaler = 1.0f / size;

        float u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        float pt_0, pt_1, pt_2 = 0.0f;
        float _attr_temp = 0.0f;
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = attr[primpoints[i]];
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        float _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0f - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) double primuv_float(
    const int geo, global double* attr, const int primnum, const double3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double _attr = 0.0;
    double size = (double)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = attr[primpoints[pt_str  ]] * (1.0 - uvw.x - uvw.y) 
              + attr[primpoints[pt_str+1]] * uvw.x
              + attr[primpoints[pt_str+2]] * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = attr[primpoints[pt_str  ]] * (1.0 - uvw.x) * (1.0 - uvw.y) 
              + attr[primpoints[pt_str+1]] * (1.0 - uvw.x) * uvw.y
              + attr[primpoints[pt_str+2]] * uvw.x         * uvw.y
              + attr[primpoints[pt_str+3]] * uvw.x         * (1.0 - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        double v_accum_scaler = 1.0 / size;

        double u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        double pt_0, pt_1, pt_2 = 0.0;
        double _attr_temp = 0.0;
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = attr[primpoints[i]];
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        double _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0 - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) float2 primuv_float2(
    const int geo, global float* attr, const int primnum, const float3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float2 _attr = (float2)(0.0f);
    float size = (float)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = vload2(primpoints[pt_str  ], attr) * (1.0f - uvw.x - uvw.y) 
              + vload2(primpoints[pt_str+1], attr) * uvw.x
              + vload2(primpoints[pt_str+2], attr) * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = vload2(primpoints[pt_str  ], attr) * (1.0f - uvw.x) * (1.0f - uvw.y) 
              + vload2(primpoints[pt_str+1], attr) * (1.0f - uvw.x) * uvw.y
              + vload2(primpoints[pt_str+2], attr) * uvw.x          * uvw.y
              + vload2(primpoints[pt_str+3], attr) * uvw.x          * (1.0f - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        float v_accum_scaler = 1.0f / size;

        float u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        float2 pt_0, pt_1, pt_2 = (float2)(0.0f);
        float2 _attr_temp = (float2)(0.0f);
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = vload2(primpoints[i], attr);
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        float _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0f - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) double2 primuv_float2(
    const int geo, global double* attr, const int primnum, const double3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double2 _attr = (double2)(0.0);
    double size = (double)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = vload2(primpoints[pt_str  ], attr) * (1.0 - uvw.x - uvw.y) 
              + vload2(primpoints[pt_str+1], attr) * uvw.x
              + vload2(primpoints[pt_str+2], attr) * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = vload2(primpoints[pt_str  ], attr) * (1.0 - uvw.x) * (1.0 - uvw.y) 
              + vload2(primpoints[pt_str+1], attr) * (1.0 - uvw.x) * uvw.y
              + vload2(primpoints[pt_str+2], attr) * uvw.x         * uvw.y
              + vload2(primpoints[pt_str+3], attr) * uvw.x         * (1.0 - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        double v_accum_scaler = 1.0 / size;

        double u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        double2 pt_0, pt_1, pt_2 = (double2)(0.0);
        double2 _attr_temp = (double2)(0.0);
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = vload2(primpoints[i], attr);
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        double _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0 - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) float3 primuv_float3(
    const int geo, global float* attr, const int primnum, const float3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float3 _attr = (float3)(0.0f);
    float size = (float)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = vload3(primpoints[pt_str  ], attr) * (1.0f - uvw.x - uvw.y) 
              + vload3(primpoints[pt_str+1], attr) * uvw.x
              + vload3(primpoints[pt_str+2], attr) * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = vload3(primpoints[pt_str  ], attr) * (1.0f - uvw.x) * (1.0f - uvw.y) 
              + vload3(primpoints[pt_str+1], attr) * (1.0f - uvw.x) * uvw.y
              + vload3(primpoints[pt_str+2], attr) * uvw.x          * uvw.y
              + vload3(primpoints[pt_str+3], attr) * uvw.x          * (1.0f - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        float v_accum_scaler = 1.0f / size;

        float u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        float3 pt_0, pt_1, pt_2 = (float3)(0.0f);
        float3 _attr_temp = (float3)(0.0f);
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = vload3(primpoints[i], attr);
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        float _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0f - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) double3 primuv_float3(
    const int geo, global double* attr, const int primnum, const double3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double3 _attr = (double3)(0.0);
    double size = (double)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = vload3(primpoints[pt_str  ], attr) * (1.0 - uvw.x - uvw.y) 
              + vload3(primpoints[pt_str+1], attr) * uvw.x
              + vload3(primpoints[pt_str+2], attr) * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = vload3(primpoints[pt_str  ], attr) * (1.0 - uvw.x) * (1.0 - uvw.y) 
              + vload3(primpoints[pt_str+1], attr) * (1.0 - uvw.x) * uvw.y
              + vload3(primpoints[pt_str+2], attr) * uvw.x         * uvw.y
              + vload3(primpoints[pt_str+3], attr) * uvw.x         * (1.0 - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        double v_accum_scaler = 1.0 / size;

        double u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        double3 pt_0, pt_1, pt_2 = (double3)(0.0);
        double3 _attr_temp = (double3)(0.0);
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = vload3(primpoints[i], attr);
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        double _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0 - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) float4 primuv_float4(
    const int geo, global float* attr, const int primnum, const float3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float4 _attr = (float4)(0.0f);
    float size = (float)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = vload4(primpoints[pt_str  ], attr) * (1.0f - uvw.x - uvw.y) 
              + vload4(primpoints[pt_str+1], attr) * uvw.x
              + vload4(primpoints[pt_str+2], attr) * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = vload4(primpoints[pt_str  ], attr) * (1.0f - uvw.x) * (1.0f - uvw.y) 
              + vload4(primpoints[pt_str+1], attr) * (1.0f - uvw.x) * uvw.y
              + vload4(primpoints[pt_str+2], attr) * uvw.x          * uvw.y
              + vload4(primpoints[pt_str+3], attr) * uvw.x          * (1.0f - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        float v_accum_scaler = 1.0f / size;

        float u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        float4 pt_0, pt_1, pt_2 = (float4)(0.0f);
        float4 _attr_temp = (float4)(0.0f);
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = vload4(primpoints[i], attr);
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        float _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0f - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) double4 primuv_float4(
    const int geo, global double* attr, const int primnum, const double3 uvw,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns an attribute value interpolated along the primitive parametric space
    using the provided uvw coordinates -> https://www.sidefx.com/docs/houdini/model/primitive_spaces.html
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - uvw: the parametric coordinates
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double4 _attr = (double4)(0.0);
    double size = (double)(pt_end - pt_str);
    if(size==3) // triangles
    {
        _attr = vload4(primpoints[pt_str  ], attr) * (1.0 - uvw.x - uvw.y) 
              + vload4(primpoints[pt_str+1], attr) * uvw.x
              + vload4(primpoints[pt_str+2], attr) * uvw.y;
    }
    else if(size==4) // quads
    {
        _attr = vload4(primpoints[pt_str  ], attr) * (1.0 - uvw.x) * (1.0 - uvw.y) 
              + vload4(primpoints[pt_str+1], attr) * (1.0 - uvw.x) * uvw.y
              + vload4(primpoints[pt_str+2], attr) * uvw.x         * uvw.y
              + vload4(primpoints[pt_str+3], attr) * uvw.x         * (1.0 - uvw.y);
    }
    else //ngons
    {
        /* for n-gons, each edge contribute 1/n to the total `u` value.
        This generates a triangle fan centerd in the polygon center.
        The strategy is find which triangle we want to interpolate by checking the
        provided `u` coords - we can then use the same triangle interpolation 
        strategy as before.
        */

        double v_accum_scaler = 1.0 / size;

        double u_scaled = floor(uvw.x * size);
        int u_scaled_start = (int)(u_scaled);
        int u_scaled_end = u_scaled_start+1;

        double4 pt_0, pt_1, pt_2 = (double4)(0.0);
        double4 _attr_temp = (double4)(0.0);
        for(int i=pt_str, j=0; i<pt_end; i++, j++)
        {
            _attr_temp = vload4(primpoints[i], attr);
            if(j==u_scaled_start) { pt_0 = _attr_temp; }
            if(j==u_scaled_end)   { pt_1 = _attr_temp; }

            pt_2 += _attr_temp * v_accum_scaler;
        }

        double _u = (uvw.x * size) - u_scaled;
        _attr = pt_0 * (1.0 - _u - uvw.y) 
              + pt_1 * _u
              + pt_2 * uvw.y;

    }
    return _attr;
}

static inline __attribute__((overloadable)) float prim_float(
    const int geo, global float* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float _attr = 0.0f;
    float size = (float)(pt_end - pt_str);

    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += attr[primpoints[i]];
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) double prim_float(
    const int geo, global double* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double _attr = 0.0;
    double size = (double)(pt_end - pt_str);

    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += attr[primpoints[i]];
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) float2 prim_float2(
    const int geo, global float* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float2 _attr = (float2)(0.0f);
    float size = (float)(pt_end - pt_str);

    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += vload2(primpoints[i], attr);
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) double2 prim_float2(
    const int geo, global double* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double2 _attr = (double2)(0.0);
    double size = (double)(pt_end - pt_str);

    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += vload2(primpoints[i], attr);
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) float3 prim_float3(
    const int geo, global float* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float3 _attr = (float3)(0.0f);
    float size = (float)(pt_end - pt_str);

    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += vload3(primpoints[i], attr);
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) double3 prim_float3(
    const int geo, global double* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double3 _attr = (double3)(0.0);
    double size = (double)(pt_end - pt_str);

    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += vload3(primpoints[i], attr);
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) float4 prim_float4(
    const int geo, global float* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    float4 _attr = (float4)(0.0f);
    float size = (float)(pt_end - pt_str);
    
    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += vload4(primpoints[i], attr);
    }

    return _attr / size;
}

static inline __attribute__((overloadable)) double4 prim_float4(
    const int geo, global double* attr, const int primnum,
    global int* primpoints_index, global int* primpoints)
{
    /* Returns a point attribute value in the center of the primitive, 
    as in averaged between all the point values. Equivalent to primuv function
    with {0.5, 0.5, 0.0} coordinates.
    
    - geo: ignored, used to maintain consistency with vex code
    - attr: point attribute array
    - primnum: the primitive index, usually `idx`
    - primpoints_index: array of primpoints indexes, used to find the actual primpoint
        values
    - primpoints: array of primitive point ids. Expected to be setup in vex beforehand
    */

    int pt_str = primpoints_index[primnum];
    int pt_end = primpoints_index[primnum+1];

    double4 _attr = (double4)(0.0);
    double size = (double)(pt_end - pt_str);
    
    for(int i=pt_str; i<pt_end; i++)
    {
        _attr += vload4(primpoints[i], attr);
    }

    return _attr / size;
}

#endif
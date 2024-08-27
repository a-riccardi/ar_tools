#ifndef __geo_quad__
#define __geo_quad__

/* Library to extract geometric information from quads. Offers overloads for both float and double precision

Methods:
    - quad_perimeter: returns the perimeter of the A, B, C, D quad
    - quad_quality: given points A, B, C, D, returns the quad quality. Optionally returns the perimeter as well
*/

static inline __attribute__((overloadable)) float quad_perimeter(const float3 A, const float3 B, const float3 C, const float3 D)
{
    /* Compute the perimeter of the ABCD quadrilateral. */
    
    float3 BA = B - A;
    float3 CB = C - B;
    float3 DC = D - C;
    float3 AD = A - D;

    return sqrt(dot(BA, BA))
         + sqrt(dot(CB, CB))
         + sqrt(dot(DC, DC))
         + sqrt(dot(AD, AD));
}

static inline __attribute__((overloadable)) double quad_perimeter(const double3 A, const double3 B, const double3 C, const double3 D)
{
    /* Compute the perimeter of the ABCD quadrilateral. */
    
    double3 BA = B - A;
    double3 CB = C - B;
    double3 DC = D - C;
    double3 AD = A - D;

    return sqrt(dot(BA, BA))
         + sqrt(dot(CB, CB))
         + sqrt(dot(DC, DC))
         + sqrt(dot(AD, AD));
}

static inline __attribute__((overloadable)) float quad_quality(const float3 A, const float3 B, const float3 C, const float3 D)
{
    /* Measure the quality of a quadrilateral as the ratio between the shortest edge and the
    longest diagonal. It also checks how much coplanar the two triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how much square the quad is.
    */

    float3 BA = B - A;
    float3 CB = C - B;
    float3 DC = D - C;
    float3 AD = A - D;

    float3 CA = C - A;
    float3 DB = D - B;

    float3 ABD_normal = normalize(cross(BA, -AD));
    float3 BCD_normal = normalize(cross(DC, -CB));

    float min_edge_squared = fmin(dot(BA, BA), fmin(dot(CB, CB), fmin(dot(DC, DC), dot(AD, AD))));
    float max_diagonal_squared = fmax(dot(CA, CA), dot(DB, DB));

    return (M_SQRT2 * sqrt(min_edge_squared) / sqrt(max_diagonal_squared))
         * (clamp((double)(dot(ABD_normal, BCD_normal)), 0.0, 1.0));
}

static inline __attribute__((overloadable)) double quad_quality(const double3 A, const double3 B, const double3 C, const double3 D)
{
    /* Measure the quality of a quadrilateral as the ratio between the shortest edge and the
    longest diagonal. It also checks how much coplanar the two triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how much square the quad is.
    */

    double3 BA = B - A;
    double3 CB = C - B;
    double3 DC = D - C;
    double3 AD = A - D;

    double3 CA = C - A;
    double3 DB = D - B;

    double3 ABD_normal = normalize(cross(BA, -AD));
    double3 BCD_normal = normalize(cross(DC, -CB));

    float min_edge_squared = fmin(dot(BA, BA), fmin(dot(CB, CB), fmin(dot(DC, DC), dot(AD, AD))));
    float max_diagonal_squared = fmax(dot(CA, CA), dot(DB, DB));

    return (M_SQRT2 * sqrt(min_edge_squared) / sqrt(max_diagonal_squared))
         * (clamp((double)(dot(ABD_normal, BCD_normal)), 0.0, 1.0));
}

static inline __attribute__((overloadable)) float quad_quality(const float3 A, const float3 B, const float3 C, const float3 D, float* perimeter)
{
    /* Measure the quality of a quadrilateral as the ratio between the shortest edge and the
    longest diagonal. It also checks how much coplanar the two triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how much square the quad is.
    */

    float3 BA = B - A;
    float3 CB = C - B;
    float3 DC = D - C;
    float3 AD = A - D;

    float ba = sqrt(dot(BA, BA));
    float cb = sqrt(dot(CB, CB));
    float dc = sqrt(dot(DC, DC));
    float ad = sqrt(dot(AD, AD));

    *perimeter = ba + cb + dc + ad;

    float3 CA = C - A;
    float3 DB = D - B;

    float3 ABD_normal = normalize(cross(BA, -AD));
    float3 BCD_normal = normalize(cross(DC, -CB));

    float min_edge = fmin(ba, fmin(cb, fmin(dc, ad)));
    float max_diagonal_squared = fmax(dot(CA, CA), dot(DB, DB));

    return (M_SQRT2 * min_edge / sqrt(max_diagonal_squared))
         * (clamp((double)(dot(ABD_normal, BCD_normal)), 0.0, 1.0));
}

static inline __attribute__((overloadable)) double quad_quality(const double3 A, const double3 B, const double3 C, const double3 D, double* perimeter)
{
    /* Measure the quality of a quadrilateral as the ratio between the shortest edge and the
    longest diagonal. It also checks how much coplanar the two triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how much square the quad is.
    */

    double3 BA = B - A;
    double3 CB = C - B;
    double3 DC = D - C;
    double3 AD = A - D;

    double ba = sqrt(dot(BA, BA));
    double cb = sqrt(dot(CB, CB));
    double dc = sqrt(dot(DC, DC));
    double ad = sqrt(dot(AD, AD));

    *perimeter = ba + cb + dc + ad;

    double3 CA = C - A;
    double3 DB = D - B;

    double3 ABD_normal = normalize(cross(BA, -AD));
    double3 BCD_normal = normalize(cross(DC, -CB));

    double min_edge = fmin(ba, fmin(cb, fmin(dc, ad)));
    double max_diagonal_squared = fmax(dot(CA, CA), dot(DB, DB));

    return (M_SQRT2 * min_edge / sqrt(max_diagonal_squared))
         * (clamp((double)(dot(ABD_normal, BCD_normal)), 0.0, 1.0));
}

#endif
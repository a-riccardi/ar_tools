#ifndef __geo_quad__
#define __geo_quad__

/* Library to extract geometric information from quads. Offers overloads for both float and double precision

Methods:
    - quad_perimeter: returns the perimeter of the A, B, C, D quad
    - quad_quality: given points A, B, C, D, returns the quad quality. Optionally returns the perimeter as well
*/

#include "math.h"

function float quad_perimeter(vector A; vector B; vector C; vector D)
{
    /* Compute the perimeter of the ABCD quadrilateral. */
    
    vector BA = B - A;
    vector CB = C - B;
    vector DC = D - C;
    vector AD = A - D;

    return sqrt(dot(BA, BA))
         + sqrt(dot(CB, CB))
         + sqrt(dot(DC, DC))
         + sqrt(dot(AD, AD));
}

function float quad_quality(vector A; vector B; vector C; vector D)
{
    /* Measure the quality of a quadrilateral as the ratio between the shortest edge and the
    longest diagonal. It also checks how much coplanar the two triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how much square the quad is.
    */

    vector BA = B - A;
    vector CB = C - B;
    vector DC = D - C;
    vector AD = A - D;

    vector CA = C - A;
    vector DB = D - B;

    vector ABD_normal = normalize(cross(BA, -AD));
    vector BCD_normal = normalize(cross(DC, -CB));
    
    float min_edge_squared = min(dot(BA, BA), min(dot(CB, CB), min(dot(DC, DC), dot(AD, AD))));
    float max_diagonal_squared = max(dot(CA, CA), dot(DB, DB));

    return (M_SQRT2 * sqrt(min_edge_squared) / sqrt(max_diagonal_squared))
         * (clamp(dot(ABD_normal, BCD_normal), 0.0, 1.0));
}

function float quad_quality(vector A; vector B; vector C; vector D; export float perimeter)
{
    /* Measure the quality of a quadrilateral as the ratio between the shortest edge and the
    longest diagonal. It also checks how much coplanar the two triangles are using the dot
    product between their normals. The result is a [0,1] value measuring how much square the quad is.
    */

    vector BA = B - A;
    vector CB = C - B;
    vector DC = D - C;
    vector AD = A - D;

    float ba = sqrt(dot(BA, BA));
    float cb = sqrt(dot(CB, CB));
    float dc = sqrt(dot(DC, DC));
    float ad = sqrt(dot(AD, AD));

    perimeter = ba + cb + dc + ad;

    vector CA = C - A;
    vector DB = D - B;

    vector ABD_normal = normalize(cross(BA, -AD));
    vector BCD_normal = normalize(cross(DC, -CB));

    float min_edge = min(ba, min(cb, min(dc, ad)));
    float max_diagonal_squared = max(dot(CA, CA), dot(DB, DB));

    return (M_SQRT2 * min_edge / sqrt(max_diagonal_squared))
         * (clamp(dot(ABD_normal, BCD_normal), 0.0, 1.0));
}

#endif
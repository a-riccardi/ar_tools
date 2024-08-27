#ifndef __geo_triangle__
#define __geo_triangle__

/* Library to extract geometric information from triangles.

Methods:
    - tris_perimeter: given points A, B, C, returns the triangle perimeter
    - tris_area: given points A, B, C, returns the triangle area
    - tris_incircle_radius: given points A, B, C, returns the triangle incircle radius
    - tris_circumcircle_radius: given points A, B, C, returns the triangle circumcircle radius
    - tris_quality: given points A, B, C, returns the triangle quality. Optionally returns the perimeter as well
*/

function float tris_perimeter(vector A; vector B; vector C)
{
    /* Returns the perimeter of the provided ABC triangle */

    vector BA = B - A;
    vector CB = C - B;
    vector AC = A - C;

    return sqrt(dot(BA, BA)) + sqrt(dot(CB, CB)) + sqrt(dot(AC, AC));
}

function float tris_area(vector A; vector B; vector C; export float abc; export float semiperimeter)
{
    /* Returns the area of the provided ABC triangle. It also returns the three sides multiplied together and the semiperimeter */

    vector BA = B - A;
    vector CB = C - B;
    vector AC = A - C;

    float ba = sqrt(dot(BA, BA));
    float cb = sqrt(dot(CB, CB));
    float ac = sqrt(dot(AC, AC));
    abc = ba * cb * ac;
    semiperimeter = (ba + cb + ac) * 0.5f;
    return sqrt(semiperimeter * (semiperimeter - ba) * (semiperimeter - cb) * (semiperimeter - ac));
}

function float tris_incircle_radius(vector A; vector B; vector C; export float area; export float semiperimeter)
{
    /* Returns the radius of the incircle of the ABC triangle. It also returns the area and semiperimeter. */

    float abc = 0.0f;
    area = tris_area(A, B, C, abc, semiperimeter);

    return area / semiperimeter;
}

function float tris_circumcircle_radius(vector A; vector B; vector C; export float area)
{
    /* Returns the radius of the incircle of the ABC triangle. It also returns the area and semiperimeter. */

    float abc = 0.0f;
    float semiperimeter = 0.0f;
    area = tris_area(A, B, C, abc, semiperimeter);

    return abc / (4.0f * area);
}

function float tris_quality(vector A; vector B; vector C)
{
    /* Measure the quality of the triangle as the ratio between incircle and circumcircle.
    The result is a [0,1] value measuring how much the triangle is equilateral
    */

    float abc = 0.0f;
    float semiperimeter = 0.0f;
    float area = tris_area(A, B, C, abc, semiperimeter);

    float incircle = area/semiperimeter;
    float circumcircle = abc / (4.0f * area);

    return (2.0f * incircle) / circumcircle;
}

function float tris_quality(vector A; vector B; vector C; export float perimeter)
{
    /* Measure the quality of the triangle as the ratio between incircle and circumcircle.
    The result is a [0,1] value measuring how much the triangle is equilateral
    */

    float abc = 0.0f;
    float semiperimeter = 0.0f;
    float area = tris_area(A, B, C, abc, semiperimeter);

    perimeter = semiperimeter * 2.0f;

    float incircle = area / semiperimeter;
    float circumcircle = abc / (4.0f * area);

    return (2.0f * incircle) / circumcircle;
}

#endif
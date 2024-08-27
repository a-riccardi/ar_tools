#ifndef __geo_tris__
#define __geo_tris__

/* Library to extract geometric information from triangles. Offers overloads for both float and double precision

Methods:
    - tris_perimeter: given points A, B, C, returns the triangle perimeter
    - tris_area: given points A, B, C, returns the triangle area
    - tris_incircle_radius: given points A, B, C, returns the triangle incircle radius
    - tris_circumcircle_radius: given points A, B, C, returns the triangle circumcircle radius
    - tris_quality: given points A, B, C, returns the triangle quality. Optionally returns the perimeter as well
*/

static inline __attribute__((overloadable)) float tris_perimeter(const float3 A, const float3 B, const float3 C)
{
    /* Returns the perimeter of the provided ABC triangle */

    float3 BA = B - A;
    float3 CB = C - B;
    float3 AC = A - C;

    return sqrt(dot(BA, BA)) + sqrt(dot(CB, CB)) + sqrt(dot(AC, AC));
}

static inline __attribute__((overloadable)) double tris_perimeter(const double3 A, const double3 B, const double3 C)
{
    /* Returns the perimeter of the provided ABC triangle */

    double3 BA = B - A;
    double3 CB = C - B;
    double3 AC = A - C;

    return sqrt(dot(BA, BA)) + sqrt(dot(CB, CB)) + sqrt(dot(AC, AC));
}

static inline __attribute__((overloadable)) float tris_area(const float3 A, const float3 B, const float3 C, float* abc, float* semiperimeter)
{
    /* Returns the area of the provided ABC triangle. It also returns the three sides multiplied together and the semiperimeter */

    float3 BA = B - A;
    float3 CB = C - B;
    float3 AC = A - C;

    float ba = sqrt(dot(BA, BA));
    float cb = sqrt(dot(CB, CB));
    float ac = sqrt(dot(AC, AC));
    *abc = ba * cb * ac;
    *semiperimeter = (ba + cb + ac) * 0.5f;
    return sqrt(*semiperimeter * (*semiperimeter - ba) * (*semiperimeter - cb) * (*semiperimeter - ac));
}

static inline __attribute__((overloadable)) double tris_area(const double3 A, const double3 B, const double3 C, double* abc, double* semiperimeter)
{
    /* Returns the area of the provided ABC triangle. It also returns the three sides multiplied together and the semiperimeter */

    double3 BA = B - A;
    double3 CB = C - B;
    double3 AC = A - C;

    double ba = sqrt(dot(BA, BA));
    double cb = sqrt(dot(CB, CB));
    double ac = sqrt(dot(AC, AC));
    *abc = ba * cb * ac;
    *semiperimeter = (ba + cb + ac) * 0.5;
    return sqrt(*semiperimeter * (*semiperimeter - ba) * (*semiperimeter - cb) * (*semiperimeter - ac));
}

static inline __attribute__((overloadable)) float tris_incircle_radius(const float3 A, const float3 B, const float3 C, float* area, float* semiperimeter)
{
    /* Returns the radius of the incircle of the ABC triangle. It also returns the area and semiperimeter. */

    float abc = 0.0f;
    *area = tris_area(A, B, C, &abc, semiperimeter);

    return *area / *semiperimeter;
}

static inline __attribute__((overloadable)) double tris_incircle_radius(const double3 A, const double3 B, const double3 C, double* area, double* semiperimeter)
{
    /* Returns the radius of the incircle of the ABC triangle. It also returns the area and semiperimeter. */

    double abc = 0.0;
    *area = tris_area(A, B, C, &abc, semiperimeter);

    return *area / *semiperimeter;
}

static inline __attribute__((overloadable)) float tris_circumcircle_radius(const float3 A, const float3 B, const float3 C, float* area)
{
    /* Returns the radius of the circumcircle of the ABC triangle. It also returns the area. */
    
    float abc = 0.0f;
    float semiperimeter = 0.0f;
    *area = tris_area(A, B, C, &abc, &semiperimeter);

    return abc / (4.0 * *area);
}

static inline __attribute__((overloadable)) double tris_circumcircle_radius(const double3 A, const double3 B, const double3 C, double* area)
{
    /* Returns the radius of the circumcircle of the ABC triangle. It also returns the area. */

    double abc = 0.0;
    double semiperimeter = 0.0;
    *area = tris_area(A, B, C, &abc, &semiperimeter);

    return abc / (4.0 * *area);
}

static inline __attribute__((overloadable)) float tris_quality(const float3 A, const float3 B, const float3 C)
{
    /* Measures the quality of the triangle as the ratio between incircle and circumcircle.
    The result is a [0,1] value measuring how much the triangle is equilateral
    */

    float abc = 0.0f;
    float semiperimeter = 0.0f;
    float area = tris_area(A, B, C, &abc, &semiperimeter);

    float incircle = area / semiperimeter;
    float circumcircle = abc / (4.0f * area);

    return (2.0f * incircle) / circumcircle;
}

static inline __attribute__((overloadable)) double tris_quality(const double3 A, const double3 B, const double3 C)
{
    /* Measures the quality of the triangle as the ratio between incircle and circumcircle.
    The result is a [0,1] value measuring how much the triangle is equilateral
    */

    double abc = 0.0;
    double semiperimeter = 0.0;
    double area = tris_area(A, B, C, &abc, &semiperimeter);

    double incircle = area / semiperimeter;
    double circumcircle = abc / (4.0 * area);

    return (2.0 * incircle) / circumcircle;
}

static inline __attribute__((overloadable)) float tris_quality(const float3 A, const float3 B, const float3 C, double* perimeter)
{
    /* Measures the quality of the triangle as the ratio between incircle and circumcircle.
    The result is a [0,1] value measuring how much the triangle is equilateral
    */

    float abc = 0.0f;
    float semiperimeter = 0.0f;
    float area = tris_area(A, B, C, &abc, &semiperimeter);

    *perimeter = semiperimeter * 2.0f;

    float incircle = area / semiperimeter;
    float circumcircle = abc / (4.0f * area);

    return (2.0f * incircle) / circumcircle;
}

static inline __attribute__((overloadable)) double tris_quality(const double3 A, const double3 B, const double3 C, double* perimeter)
{
    /* Measures the quality of the triangle as the ratio between incircle and circumcircle.
    The result is a [0,1] value measuring how much the triangle is equilateral
    */

    double abc = 0.0;
    double semiperimeter = 0.0;
    double area = tris_area(A, B, C, &abc, &semiperimeter);

    *perimeter = semiperimeter * 2.0;

    double incircle = area / semiperimeter;
    double circumcircle = abc / (4.0 * area);

    return (2.0 * incircle) / circumcircle;
}

#endif
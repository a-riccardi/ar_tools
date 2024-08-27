#ifndef __geo_tracing__
#define __geo_tracing__

#define EPSILON 0.0001f
#define EPSILON_D 0.0001

/* Library to handle geometry tracing and barycentric coordinates intersection.

Methods:
    - project_on_plane: projects a point P on the plane defined by point C and normal N
    - project_on_plane: projects a point P on the plane defined by plane equation abcd (abc is the normal, d is the signed distance from origin)
    - inbound_point_triangle: checks if a point P is inside the triangle defined by A, B, C. Returns the parametric coordinates if true
*/ 

static inline __attribute__((overloadable)) float3 project_on_plane(const float3 P, const float3 C, const float3 N)
{
    /* Project point P on the plane defined by point C and normal N */

    return P - (N * dot(N, P - C));
}

static inline __attribute__((overloadable)) double3 project_on_plane(const double3 P, const double3 C, const double3 N)
{
    /* Project point P on the plane defined by point C and normal N */

    return P - (N * dot(N, P - C));
}

static inline __attribute__((overloadable)) float3 project_on_plane(const float3 P, const float4 abcd)
{
    /* Project point P on the plane defined by plane equation abcd.
    abc is expected to be the plane normal, and d is the signed distance from origin
    */

    return P - (dot(abcd.xyz, P) + abcd.w) * abcd.xyz;
}

static inline __attribute__((overloadable)) double3 project_on_plane(const double3 P, const double4 abcd)
{
    /* Project point P on the plane defined by plane equation abcd.
    abc is expected to be the plane normal, and d is the signed distance from origin
    */

    return P - (dot(abcd.xyz, P) + abcd.w) * abcd.xyz;
}

static __attribute__((overloadable)) int inbound_point_triangle(const float3 P,
    const float3 A, const float3 B, const float3 C, float3* uvw)
{
    /* Check if the provided point P falls inside the ABC triangle.       B
    The triangle is expected to be in counter-cloxkwise fashion,         / \
    with A in the lower right corner as illustrated here.               /   \
    The function also returns the barycentric uvw coordinates of P     C-----A
    
    Note: the point is expected to be coplanar with ABC - use function
    `project_on_plane` above to obtain P projection on the triangle plane
    */
    
    // vector perpendicular to triangle's plane describing the position of P relative to each edge
    float3 area_vector;
    float3 v0v1 = B - A;
    float3 v0v2 = C - A;
    float3 N = cross(v0v1, v0v2);
    
    float denom = 1.0f / max(dot(N, N), EPSILON);
    float w = 0.0f;
    float u = 0.0f;

    float3 AC_edge = A - C; 
    float3 PC = P - C;
    area_vector = cross(AC_edge, PC);
    if (dot(N, area_vector) < EPSILON)      { return 0; } // P outside of AC

    float3 BA_edge = B - A; 
    float3 PA = P - A;
    area_vector = cross(BA_edge, PA);
    if ((w = dot(N, area_vector)) < EPSILON) { return 0; } // P outside of BA
 
    float3 CB_edge = C - B; 
    float3 PB = P - B;
    area_vector = cross(CB_edge, PB);
    if ((u = dot(N, area_vector)) < EPSILON) { return 0; } // P outside of CB

    w *= denom;
    u *= denom;

    *uvw = (float3)(u, 1.0f - w - u, w);
    return 1;
}

static __attribute__((overloadable)) int inbound_point_triangle(const double3 P,
    const double3 A, const double3 B, const double3 C, double3* uvw)
{
    /* Check if the provided point P falls inside the ABC triangle.       B
    The triangle is expected to be in counter-cloxkwise fashion,         / \
    with A in the lower right corner as illustrated here.               /   \
    The function also returns the barycentric uvw coordinates of P     C-----A
    
    Note: the point is expected to be coplanar with ABC - use function
    `project_on_plane` above to obtain P projection on the triangle plane
    */
    
    // vector perpendicular to triangle's plane describing the position of P relative to each edge
    double3 area_vector;
    double3 v0v1 = B - A;
    double3 v0v2 = C - A;
    double3 N = cross(v0v1, v0v2);
    
    double denom = 1.0f / max(dot(N, N), EPSILON_D);
    double w = 0.0f;
    double u = 0.0f;

    double3 AC_edge = A - C; 
    double3 PC = P - C;
    area_vector = cross(AC_edge, PC);
    if (dot(N, area_vector) < EPSILON_D)      { return 0; } // P outside of AC

    double3 BA_edge = B - A; 
    double3 PA = P - A;
    area_vector = cross(BA_edge, PA);
    if ((w = dot(N, area_vector)) < EPSILON_D) { return 0; } // P outside of BA
 
    double3 CB_edge = C - B; 
    double3 PB = P - B;
    area_vector = cross(CB_edge, PB);
    if ((u = dot(N, area_vector)) < EPSILON_D) { return 0; } // P outside of CB

    w *= denom;
    u *= denom;

    *uvw = (double3)(u, 1.0f - w - u, w);
    return 1;
}

#endif
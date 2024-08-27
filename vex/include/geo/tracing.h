#ifndef __geo_tracing__
#define __geo_tracing__

/* Library to handle geometry tracing and barycentric coordinates intersection.

Methods:
    - project_on_plane: projects a point P on the plane defined by point C and normal N
    - project_on_plane: projects a point P on the plane defined by plane equation abcd (abc is the normal, d is the signed distance from origin)
    - inbound_point_triangle: checks if a point P is inside the triangle defined by A, B, C. Returns the parametric coordinates if true
*/

#define EPSILON 0.0001f

function vector project_on_plane(const vector P; const vector C; const vector N)
{
    /* Project point P on the plane defined by point C and normal N */

    return P - (N * dot(N, P - C));
}

function vector project_on_plane(const vector P; const vector4 abcd)
{
    /* Project point P on the plane defined by plane equation abcd.
    abc is expected to be the plane normal, and d is the signed distance from origin
    */

    return P - (dot(vector(abcd), P) + abcd.w) * vector(abcd);
}

function int inbound_point_triangle(const vector P;
    const vector A; const vector B; const vector C; export vector uvw)
{
    /* Check if the provided point P falls inside the ABC triangle.       B
    The triangle is expected to be in counter-clockwise fashion,         / \
    with A in the lower right corner as illustrated here.               /   \
    The function also returns the barycentric uvw coordinates of P     C-----A
    
    Note: the point is expected to be coplanar with ABC - use function
    `project_on_plane` above to obtain P projection on the triangle plane
    */
    
    // vector perpendicular to triangle's plane describing the position of P relative to each edge
    vector area_vector;
    vector v0v1 = B - A;
    vector v0v2 = C - A;
    vector N = cross(v0v1, v0v2);
    
    float denom = 1.0f / max(dot(N, N), EPSILON);
    float w = 0.0f;
    float u = 0.0f;

    vector AC_edge = A - C; 
    vector PC = P - C;
    area_vector = cross(AC_edge, PC);
    if (dot(N, area_vector) < EPSILON)      { return 0; } // P outside of AC

    vector BA_edge = B - A; 
    vector PA = P - A;
    area_vector = cross(BA_edge, PA);
    if ((w = dot(N, area_vector)) < EPSILON) { return 0; } // P outside of BA
 
    vector CB_edge = C - B; 
    vector PB = P - B;
    area_vector = cross(CB_edge, PB);
    if ((u = dot(N, area_vector)) < EPSILON) { return 0; } // P outside of CB

    w *= denom;
    u *= denom;

    uvw = set(u, 1.0f - w - u, w);
    return 1;
}

#endif
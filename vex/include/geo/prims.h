#ifndef __geo_prims__
#define __geo_prims__

#ifndef PRIM_UVW_CENTER
    #define PRIM_UVW_CENTER {0.5f, 0.5f, 0.0f}
#endif

/*  Library containing primitive-related utilities, such as gather points in an ordered fashion or collect neighbours

Methods:
    - neighbours_prim: given a primitive number, returns an array of neighbouring primitive indexes
    - primpoints_ordered: given a primitive number, returns a list of the primitive points as they appear following the
                          primitive edge, starting with point 0.
*/

function int[] neighbours_prim(int geo; int primnum)
{
    /* Returns a list of primitive that share an edge with the given primnum */

    int prim_neighbours[];
    foreach(int vtx; primvertices(geo, primnum))
    {
        int this_hedge = vertexhedge(geo, vtx);
        
        if (hedge_isvalid(geo, this_hedge))
        {
            int equiv_hedge = this_hedge;
            
            do
            {
                equiv_hedge = hedge_nextequiv(geo, equiv_hedge);
                
                if (equiv_hedge != this_hedge)
                {
                    push(prim_neighbours, hedge_prim(geo, equiv_hedge));   
                }
            
            } while (equiv_hedge != this_hedge);
        }
    }
    
    return prim_neighbours;
}

function int[] primpoints_ordered(int geo; int primnum; int primpoints[])
{
    /* Returns an array of point indexes for the given primitive, ordered along the prim edges.
    Consider the two quads on the right, it's easy to       A---B       A---B
    see how can be complicated to build proper normals      |   |       |   |
    and extract shape info if the points are not ordered    D---C       C---D
    in a consistent manner. This function starts at A and returns all the point in the order they appear
    along the edges, allowing for consistent code approach. */

    int pt = primpoints[0];
    int h_edge = pointhedge(geo, pt);
    int primpoints_ordered[];
    resize(primpoints_ordered, len(primpoints));

    // ensure that the half-edge that we got belongs to this primitive
    while(hedge_prim(geo, h_edge) != primnum)
    {
        h_edge = pointhedgenext(geo, h_edge);
    }

    int i = 0;
    do // until we cycled through all edges
    {
        primpoints_ordered[i] = pt;
        
        h_edge = hedge_next(geo, h_edge);
        pt = hedge_srcpoint(geo, h_edge);
        i++;
    }
    while(pt != primpoints[0]);

    return primpoints_ordered;
}

#endif

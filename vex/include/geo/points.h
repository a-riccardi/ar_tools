#ifndef __geo_points__
#define __geo_points__

/*  Library containing point-related utilities, such as check unshared points

Methods:
    - is_unshared: given a point nnumber, returns 1 if the point is unshared, 0 otherwise
*/

function int is_unshared(int geo; int ptnum)
{
    /* Returns 1 if the point is part of an unshared edge, 0 otherwise */

    return neighbourcount(geo, ptnum) != len(pointprims(geo, ptnum));
}

#endif
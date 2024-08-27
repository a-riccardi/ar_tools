/* Utility library to work with arrays */

#ifndef __types_array__
#define __types_array__

function void intersect(string a[]; string b[]; export string result[])
{
    /* Computes the intersection between string arrays a and b and return the result */
    
    int a_over_b = len(a) <= len(b);
    result =  a_over_b ? a : b;

    if(a_over_b)
    {
        foreach(string _a; result)
        {
            int index = find(b, _a);
            if(index < 0) { removevalue(result, _a); }
        }
    }
    else
    {
        foreach(string _b; result)
        {
            int index = find(a, _b);
            if(index < 0) { removevalue(result, _b); }
        }
    }
}

function void shuffle(int a[]; float seed)
{
    /* In-place shuffle array components in a random fashion */

    int size = len(a);
    int src_pointer, dst_pointer, swap;
    vector2 random_pointers;
    for(int i=0; i<size; i++)
    {
        // generate a couple of randomly selected array indexes to swap
        random_pointers = rand((i * 3.1415f) + seed);
        random_pointers = floor(random_pointers * size);

        // cast to int
        src_pointer = (int)random_pointers.y; 
        dst_pointer = (int)random_pointers.x; 

        // perform swap
        swap = a[dst_pointer];
        a[dst_pointer] = a[src_pointer];
        a[src_pointer] = swap;
    }
}

function void merge(export int dst[]; int src[])
{
    /* Merge src into dst ensuring no duplicate items are added to dst.
    The function is optimized for small src incrementally added to a growing dst array
    */

    int add[];
    resize(add, len(dst), 1);
    foreach(int dst_item; dst)
    {
        for(int i=0; i<len(src); i++)
        {
            if(add[i] && src[i] == dst_item) { add[i] = 0; }
        }
    }

    for(int i=0; i<len(src); i++)
    {
        if(add[i]) { append(dst, src[i]); }
    }
}

#endif
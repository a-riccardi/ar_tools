/* This module contains general utility functions. The transfer_attrib functions help transfer attributes
from a geometry piece to another, regardless of attribute class or type
*/

#ifndef __ar_utils__
#define __ar_utils__

#define FEQUAL_EPSILON 0.00001

function float   greater(const float   a; const float   b) { return a > b; }
function vector2 greater(const vector2 a; const vector2 b) { return set(a.x > b.x, a.y > b.y); }
function vector  greater(const vector  a; const vector  b) { return set(a.x > b.x, a.y > b.y, a.z > b.z); }
function vector4 greater(const vector4 a; const vector4 b) { return set(a.x > b.x, a.y > b.y, a.z > b.z, a.w > b.w); }


function float   greater_equal(const float   a; const float   b) { return a >= b; }
function vector2 greater_equal(const vector2 a; const vector2 b) { return set(a.x >= b.x, a.y >= b.y); }
function vector  greater_equal(const vector  a; const vector  b) { return set(a.x >= b.x, a.y >= b.y, a.z >= b.z); }
function vector4 greater_equal(const vector4 a; const vector4 b) { return set(a.x >= b.x, a.y >= b.y, a.z >= b.z, a.w >= b.w); }
/*
function float sum(float   v) { return v; }
function float sum(vector2 v) { return v.x + v.y; }
function float sum(vector  v) { return v.x + v.y + v.z; }
function float sum(vector4 v) { return v.x + v.y + v.z + v.w; }
*/
function int any(const float   v) { return int(greater(sum(v), 0.0f)); }
function int any(const vector2 v) { return int(greater(sum(v), 0.0f)); }
function int any(const vector  v) { return int(greater(sum(v), 0.0f)); }
function int any(const vector4 v) { return int(greater(sum(v), 0.0f)); }


int fequal(float a; float b) { return abs(a-b) < FEQUAL_EPSILON; }

void transfer_attrib(int src_geo; string src_class; int src_id;
    int dst_geo; string dst_class; int dst_id; string attribname;
    int attribsize; int attribtype; export dict out; int silent_fail)
{
    /* Transfer an attribute from a source geometry, class and id to an arbitrary destination geometry, class and id.
    Exports attribsize, attribtype and attribvalue in "S", "T" and "value" keys of 'out' dict

    if 'silent_fail' is > 0, ignores non-existing attributes, raises an error otherwise
    */

    out["S"] = attribsize;
    out["T"] = attribtype;
    if(attribtype<0 || attribsize<1)
    { 
        if(silent_fail > 0) { return; }
        error("Invalid Type[%i] and/or Size[%i] for attrib '%s'", attribtype, attribsize, attribname);
    }

    if(attribtype==0) // Integer
    {
        int v = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype;
    }
    else if (attribtype==1) // Float/Vector
    {
        if(attribsize==1) // Float
        {
            float v = attrib(src_geo, src_class, attribname, src_id);
            setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
            out["value"] = v;
            out["S"] = attribsize;
            out["T"] = attribtype;
        }
        else if(attribsize==2) // Vector2
        {
            vector2 v = attrib(src_geo, src_class, attribname, src_id);
            setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
            out["value"] = v;
            out["S"] = attribsize;
            out["T"] = attribtype;
        }
        else if(attribsize==3) // Vector
        {
            vector v = attrib(src_geo, src_class, attribname, src_id);
            setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
            out["value"] = v;
            out["S"] = attribsize;
            out["T"] = attribtype;
        }
        else if(attribsize==4) // Vector4 or 2x2 Matrix
        {
            vector4 v = attrib(src_geo, src_class, attribname, src_id);
            setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
            out["value"] = v;
            out["S"] = attribsize;
            out["T"] = attribtype;
        }
        else if(attribsize==9) // 3x3 Matrix
        {
            matrix3 v = attrib(src_geo, src_class, attribname, src_id);
            setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
            out["value"] = v;
            out["S"] = attribsize;
            out["T"] = attribtype;
        }
        else if(attribsize==16) // 4x4 Matrix
        {
            matrix v = attrib(src_geo, src_class, attribname, src_id);
            setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
            out["value"] = v;
            out["S"] = attribsize;
            out["T"] = attribtype;
        }
    }
    else if (attribtype==2) // String
    {
        string v = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype;
    }
    else if (attribtype==3) // Array of Integers
    {
        int v[] = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype;
    }
    else if (attribtype==4) // Array of Floats
    {
        float v[] = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype;
    }
    else if (attribtype==5) // Array of Strings
    {
        string v[] = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype; 
    }
    else if (attribtype==6) // Dictionary
    {
        dict v = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype;
    }
    else if (attribtype==7) // Array of Dictionaries
    {
        dict v[] = attrib(src_geo, src_class, attribname, src_id);
        setattrib(dst_geo, dst_class, attribname, dst_id, 0, v);
        out["value"] = v;
        out["S"] = attribsize;
        out["T"] = attribtype;
    }
}

void transfer_attrib(int src_geo; string src_class; int src_id;
    int dst_geo; string dst_class; int dst_id; string attribname)
{
    dict out;
    int S = attribsize(src_geo, src_class, attribname);
    int T = attribtype(src_geo, src_class, attribname);

    transfer_attrib(src_geo, src_class, src_id,
        dst_geo, dst_class, dst_id, attribname,
        S, T, out, 1);
}

void transfer_attrib(int src_geo; string src_class; int src_id;
    int dst_geo; string dst_class; int dst_id; string attribname;
    export dict out)
{
    int S = attribsize(src_geo, src_class, attribname);
    int T = attribtype(src_geo, src_class, attribname);

    transfer_attrib(src_geo, src_class, src_id,
        dst_geo, dst_class, dst_id, attribname,
        S, T, out, 1);
}

void transfer_attrib_strict(int src_geo; string src_class; int src_id;
    int dst_geo; string dst_class; int dst_id; string attribname)
{
    dict out;
    int S = attribsize(src_geo, src_class, attribname);
    int T = attribtype(src_geo, src_class, attribname);

    transfer_attrib(src_geo, src_class, src_id,
        dst_geo, dst_class, dst_id, attribname,
        S, T, out, 0);
}

void transfer_attrib_strict(int src_geo; string src_class; int src_id;
    int dst_geo; string dst_class; int dst_id; string attribname;
    export dict out)
{
    int S = attribsize(src_geo, src_class, attribname);
    int T = attribtype(src_geo, src_class, attribname);

    transfer_attrib(src_geo, src_class, src_id,
        dst_geo, dst_class, dst_id, attribname,
        S, T, out, 0);
}

#endif
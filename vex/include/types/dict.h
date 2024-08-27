/* Utility library to work with dictionaries */

#ifndef __types_dict__
#define __types_dict__

function int pop(dict d; string key)
{
    int _v = d[key];
    removeindex(d, key);
    return _v;
}

function float pop(dict d; string key)
{
    float _v = d[key];
    removeindex(d, key);
    return _v;
}

function vector2 pop(dict d; string key)
{
    vector2 _v = d[key];
    removeindex(d, key);
    return _v;
}

function vector pop(dict d; string key)
{
    vector _v = d[key];
    removeindex(d, key);
    return _v;
}

function vector4 pop(dict d; string key)
{
    vector4 _v = d[key];
    removeindex(d, key);
    return _v;
}

function string pop(dict d; string key)
{
    string _v = d[key];
    removeindex(d, key);
    return _v;
}

function int[] pop(dict d; string key)
{
    int _v[] = d[key];
    removeindex(d, key);
    return _v;
}

function float[] pop(dict d; string key)
{
    float _v[] = d[key];
    removeindex(d, key);
    return _v;
}

function string[] pop(dict d; string key)
{
    string _v[] = d[key];
    removeindex(d, key);
    return _v;
}

function dict pop(dict d; string key)
{
    dict _v = d[key];
    removeindex(d, key);
    return _v;
}

function dict[] pop(dict d; string key)
{
    dict _v[] = d[key];
    removeindex(d, key);
    return _v;
}

#endif
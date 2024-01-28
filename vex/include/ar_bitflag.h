/* Utilities to work with bitmasks and bitshift values */

#ifndef __ar_bitflag__
#define __ar_bitflag__

function void set_flag(export int bitmask; const int flag )
{
    bitmask |= flag;
}

function void remove_flag(export int bitmask; const int flag)
{
    bitmask &= ~flag;
}

function int is_flag_set(int bitmask; const int flag)
{
    return (bitmask & flag) == flag;
}

#endif
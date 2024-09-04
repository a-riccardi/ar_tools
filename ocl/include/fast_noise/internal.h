#ifndef __fast_noise_internal__
#define __fast_noise_internal__

/* This module contains the actual implementation of the noise called in fast_noise/api.h
The general module structure is as follows:
    -definition of utils functions
    -definition of basic noise func
    -macro for fractal generation
    -macro for cellular function
    -actual declaration with macro usage
*/

#define FAST_TO_INT(x) *(int*)(&##x)

__constant float HASH_TO_FLOAT = 1.0f / 2147483648.0f;
__constant int X_PRIME = 1619;
__constant int Y_PRIME = 31337;
__constant int Z_PRIME = 6791;
__constant float3 PRIME_VECTOR = {1619, 31337, 6791};
__constant int BIT10_MASK = 1023;
__constant float N_511_5 = 511.5f;
__constant float F3 = 1.0f / 3.0f;
__constant float G3 = 1.0f / 6.0f;
__constant float G33 = (3.0f / 6.0f) - 1.0f;
__constant float N_0_6 = 0.6f;
__constant float3 F3_V = {1.0f/3.0f, 1.0f/3.0f, 1.0f/3.0f};
__constant float3 G3_V = {1.0f/6.0f, 1.0f/6.0f, 1.0f/6.0f};
__constant float3 G33_V = {(3.0f / 6.0f) - 1.0f, (3.0f / 6.0f) - 1.0f, (3.0f / 6.0f) - 1.0f};
__constant float CUBIC_BOUNDING = 1.0f / (1.5f * 1.5f * 1.5f);

__constant float X_GRAD[] = { 0,-1, 0, 1, 0, 0, 0, 0,-1, 1,-1, 1,-1, 1,-1, 1 }; 
__constant float Y_GRAD[] = {-1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0,-1,-1, 1, 1 };
__constant float Z_GRAD[] = {-1, 0, 1, 0,-1,-1, 1, 1,-1,-1, 1, 1, 0, 0, 0, 0 };

/* add all components of a vector together and return the result */
static inline float sum (const float  a) { return a; }
static inline float sum2(const float2 a) { return a.x + a.y; }
static inline float sum3(const float3 a) { return a.x + a.y + a.z; }
static inline float sum4(const float4 a) { return a.x + a.y + a.z + a.w; }

static inline float lerp_cubic(const float a, const float b, const float c, const float d, const float t)
{
    float p = (d-c) - (a-b);
    return (t * t * t * p) + (t * t * (a - b - p)) + (t * (c - a) + b);
}

static inline float interp_quintic1(const float t)
{
    return  t * t * t *((t * 6.0f - 15.0f) * t + 10.0f);
}

static inline float3 interp_quintic3(const float3 t)
{
    return  t * t * t *((t * 6.0f - 15.0f) * t + 10.0f);
}

static inline int hash_vector(const float3 P, const int seed)
{   
    float3 _P = floor(P);
    int hash = seed;

    hash ^= (int)(_P.x);
    hash ^= (int)(_P.y);
    hash ^= (int)(_P.z);

    hash = hash * hash * 60493 * hash;

    return (hash >> 13) ^ hash;
}

static inline int hash_coords(const float x, const float y, const float z, const int seed)
{
    int hash = seed;

    hash ^= (int)rint(x);
    hash ^= (int)rint(y);
    hash ^= (int)rint(z);

    hash = hash * hash * 60493 * hash;

    return (hash >> 13) ^ hash;
}

static inline int hash_vector_hb(const float3 P, const int seed)
{
    float3 _P = rint(P);
    int hash = seed;

    hash ^= (int)(_P.x);
    hash ^= (int)(_P.y);
    hash ^= (int)(_P.z);

    return hash * hash * 60493 * hash;
}

static inline int hash_coords_hb(const float x, const float y, const float z, const int seed)
{
    int hash = seed;

    hash ^= (int)rint(x);
    hash ^= (int)rint(y);
    hash ^= (int)rint(z);

    return hash * hash * 60493 * hash;
}

static inline float value_coord_vector(const float3 P, const int seed)
{
    float3 _P = rint(P);
    int hash = seed;

    hash ^= (int)(_P.x);
    hash ^= (int)(_P.y);
    hash ^= (int)(_P.z);

    hash = hash * hash * 60493 * hash;

    return (float)(hash) * HASH_TO_FLOAT;
}

static inline float value_coord_coords(const float x, const float y, const float z, const int seed)
{
    int hash = seed;

    hash ^= (int)rint(x);
    hash ^= (int)rint(y);
    hash ^= (int)rint(z);

    hash = hash * hash * 60493 * hash;

    return (float)(hash) * HASH_TO_FLOAT;
}

static inline float gradient_coord_coords(const float i, const float j, const float k, const float x, const float y, const float z, const int seed)
{
    int hash = hash_coords(i, j, k, seed);
    int idx = hash & 7;

    return x * (float)(X_GRAD[idx]) + y * (float)(Y_GRAD[idx]) + z * (float)(Z_GRAD[idx]);
}

static inline float gradient_coord_vector(const float3 ijk, const float3 P, const int seed)
{
    int hash = hash_vector(ijk, seed);
    int idx = hash & 7;

    return P.x * (float)(X_GRAD[idx]) + P.y * (float)(Y_GRAD[idx]) + P.z * (float)(Z_GRAD[idx]);
}

static inline float white_noise_single(const float3 P, const int seed)
{
    int x = (int)rint(P.x);
    int y = (int)rint(P.y);
    int z = (int)rint(P.z);
    
    return value_coord_coords(
        (x ^ (x >> 16)) * X_PRIME,
        (y ^ (y >> 16)) * Y_PRIME,
        (z ^ (z >> 16)) * Z_PRIME,
        seed
    );
}

static inline float value_single(const float3 P, const int seed)
{
    float3 P_floor = floor(P);
    float3 P0 = P_floor * PRIME_VECTOR;
    float3 P1 = P0 + PRIME_VECTOR;
    float3 Ps = interp_quintic3(P - P_floor);
    
    return mix(
        mix(
            mix(value_coord_coords(P0.x, P0.y, P0.z, seed), value_coord_coords(P1.x, P0.y, P0.z, seed), Ps.x),
            mix(value_coord_coords(P0.x, P1.y, P0.z, seed), value_coord_coords(P1.x, P1.y, P0.z, seed), Ps.x),
            Ps.y
        ),
        mix(
            mix(value_coord_coords(P0.x, P0.y, P1.z, seed), value_coord_coords(P1.x, P0.y, P1.z, seed), Ps.x),
            mix(value_coord_coords(P0.x, P1.y, P1.z, seed), value_coord_coords(P1.x, P1.y, P1.z, seed), Ps.x),
            Ps.y
        ),
        Ps.z
    );
}

static inline float perlin_single(const float3 P, const int seed)
{
    float3 P_floor = floor(P);
    float3 P0 = P_floor * PRIME_VECTOR;
    float3 P1 = P0 + PRIME_VECTOR;

    float3 f0 = P - P_floor;
    float3 f1 = f0 - 1.0f; 

    P_floor = f0;
    P_floor = interp_quintic3(P_floor);

    return mix(
		mix(
			mix(gradient_coord_coords(P0.x, P0.y, P0.z, f0.x, f0.y, f0.z, seed), gradient_coord_coords(P1.x, P0.y, P0.z, f1.x, f0.y, f0.z, seed), P_floor.x),
			mix(gradient_coord_coords(P0.x, P1.y, P0.z, f0.x, f1.y, f0.z, seed), gradient_coord_coords(P1.x, P1.y, P0.z, f1.x, f1.y, f0.z, seed), P_floor.x), P_floor.y),
		mix(
			mix(gradient_coord_coords(P0.x, P0.y, P1.z, f0.x, f0.y, f1.z, seed), gradient_coord_coords(P1.x, P0.y, P1.z, f1.x, f0.y, f1.z, seed), P_floor.x),
			mix(gradient_coord_coords(P0.x, P1.y, P1.z, f0.x, f1.y, f1.z, seed), gradient_coord_coords(P1.x, P1.y, P1.z, f1.x, f1.y, f1.z, seed), P_floor.x), P_floor.y), P_floor.z);
}

static inline float simplex_single(const float3 P, const int seed)
{
    float f = sum3(P) * F3;
    float3 P0 = floor(P + f); 

    float3 ijk = P0 * PRIME_VECTOR;

    float g = sum3(P0) * G3;
    P0 = P - (P0 - g);

    int x0_ge_y0 = P0.x >= P0.y;
	int y0_ge_z0 = P0.y >= P0.z;
	int x0_ge_z0 = P0.x >= P0.z; 

    float3 ijk1 = (float3)(
        x0_ge_y0 & x0_ge_z0,
        (!(x0_ge_y0)) & y0_ge_z0,
        (!(x0_ge_z0)) & !(y0_ge_z0)
    );

	float3 ijk2 = (float3)(
        x0_ge_y0 | x0_ge_z0,
        (!(x0_ge_y0)) | y0_ge_z0,
        !(x0_ge_z0 & y0_ge_z0)
    );

    float3 xyz1 = G3_V + (P0 - ijk1);
    float3 xyz2 = F3_V + (P0 - ijk2);
	float3 xyz3 = P0 + G33_V;
	
    float4 x4 = (float4)(P0.x, xyz1.x, xyz2.x, xyz3.x);
    float4 y4 = (float4)(P0.y, xyz1.y, xyz2.y, xyz3.y);
    float4 z4 = (float4)(P0.z, xyz1.z, xyz2.z, xyz3.z);
    float4 t =  (float4)(0.6f, 0.6f, 0.6f, 0.6f) - (x4 * x4) - (y4 * y4) - (z4 * z4);

    int4 n = t >= 0.0f;

    t *= t * t * t;

    float4 v = t * (float4)(
        gradient_coord_vector(ijk, P0, seed),
        gradient_coord_vector(ijk + (PRIME_VECTOR * ijk1), xyz1, seed),
        gradient_coord_vector(ijk + (PRIME_VECTOR * ijk2), xyz2, seed),
        gradient_coord_vector(ijk + PRIME_VECTOR, xyz3, seed)
    );

	return sum4(v * as_float4(n)) * 32.0f;
}

static inline float cubic_single(const float3 P, const int seed)
{
    float3 P_floor = floor(P);

    float3 P1 = P_floor * PRIME_VECTOR;
    float3 P0 = P1 - PRIME_VECTOR;
    float3 P2 = P1 + PRIME_VECTOR;
    float3 P3 = P2 + PRIME_VECTOR;

    float3 Ps = P - P_floor;

	return CUBIC_BOUNDING * lerp_cubic(
		lerp_cubic(
			lerp_cubic(value_coord_coords(P0.x, P0.y, P0.z, seed), value_coord_coords(P1.x, P0.y, P0.z, seed), value_coord_coords(P2.x, P0.y, P0.z, seed), value_coord_coords(P3.x, P0.y, P0.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P1.y, P0.z, seed), value_coord_coords(P1.x, P1.y, P0.z, seed), value_coord_coords(P2.x, P1.y, P0.z, seed), value_coord_coords(P3.x, P1.y, P0.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P2.y, P0.z, seed), value_coord_coords(P1.x, P2.y, P0.z, seed), value_coord_coords(P2.x, P2.y, P0.z, seed), value_coord_coords(P3.x, P2.y, P0.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P3.y, P0.z, seed), value_coord_coords(P1.x, P3.y, P0.z, seed), value_coord_coords(P2.x, P3.y, P0.z, seed), value_coord_coords(P3.x, P3.y, P0.z, seed), Ps.x),
			Ps.y),
		lerp_cubic(
			lerp_cubic(value_coord_coords(P0.x, P0.y, P1.z, seed), value_coord_coords(P1.x, P0.y, P1.z, seed), value_coord_coords(P2.x, P0.y, P1.z, seed), value_coord_coords(P3.x, P0.y, P1.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P1.y, P1.z, seed), value_coord_coords(P1.x, P1.y, P1.z, seed), value_coord_coords(P2.x, P1.y, P1.z, seed), value_coord_coords(P3.x, P1.y, P1.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P2.y, P1.z, seed), value_coord_coords(P1.x, P2.y, P1.z, seed), value_coord_coords(P2.x, P2.y, P1.z, seed), value_coord_coords(P3.x, P2.y, P1.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P3.y, P1.z, seed), value_coord_coords(P1.x, P3.y, P1.z, seed), value_coord_coords(P2.x, P3.y, P1.z, seed), value_coord_coords(P3.x, P3.y, P1.z, seed), Ps.x),
			Ps.y),
		lerp_cubic(
			lerp_cubic(value_coord_coords(P0.x, P0.y, P2.z, seed), value_coord_coords(P1.x, P0.y, P2.z, seed), value_coord_coords(P2.x, P0.y, P2.z, seed), value_coord_coords(P3.x, P0.y, P2.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P1.y, P2.z, seed), value_coord_coords(P1.x, P1.y, P2.z, seed), value_coord_coords(P2.x, P1.y, P2.z, seed), value_coord_coords(P3.x, P1.y, P2.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P2.y, P2.z, seed), value_coord_coords(P1.x, P2.y, P2.z, seed), value_coord_coords(P2.x, P2.y, P2.z, seed), value_coord_coords(P3.x, P2.y, P2.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P3.y, P2.z, seed), value_coord_coords(P1.x, P3.y, P2.z, seed), value_coord_coords(P2.x, P3.y, P2.z, seed), value_coord_coords(P3.x, P3.y, P2.z, seed), Ps.x),
			Ps.y),
		lerp_cubic(
			lerp_cubic(value_coord_coords(P0.x, P0.y, P3.z, seed), value_coord_coords(P1.x, P0.y, P3.z, seed), value_coord_coords(P2.x, P0.y, P3.z, seed), value_coord_coords(P3.x, P0.y, P3.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P1.y, P3.z, seed), value_coord_coords(P1.x, P1.y, P3.z, seed), value_coord_coords(P2.x, P1.y, P3.z, seed), value_coord_coords(P3.x, P1.y, P3.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P2.y, P3.z, seed), value_coord_coords(P1.x, P2.y, P3.z, seed), value_coord_coords(P2.x, P2.y, P3.z, seed), value_coord_coords(P3.x, P2.y, P3.z, seed), Ps.x),
			lerp_cubic(value_coord_coords(P0.x, P3.y, P3.z, seed), value_coord_coords(P1.x, P3.y, P3.z, seed), value_coord_coords(P2.x, P3.y, P3.z, seed), value_coord_coords(P3.x, P3.y, P3.z, seed), Ps.x),
			Ps.y),
		Ps.z);
}

/* Gradient perturbation */

#define GRADIENT_COORD(_x,_y,_z) \
    int hash##_x##_y##_z = hash_coords_hb(P##_x.x, P##_y.y, P##_z.z, seed);\
    float3 P##_x##_y##_z = (float3)(\
        (float)(hash##_x##_y##_z & BIT10_MASK),\
        (float)((hash##_x##_y##_z >> 10) & BIT10_MASK),\
        (float)((hash##_x##_y##_z >> 20) & BIT10_MASK) \
    ) - N_511_5;

static inline float3 perturb_gradient_single(float3 P, int seed,  float p_amplitude, float p_frequency)
{
    float3 P_scaled = P * p_frequency;
    float3 P_floor = floor(P_scaled);
    float3 P0 = P_floor * PRIME_VECTOR;
    float3 P1 = P0 + PRIME_VECTOR;

    P_floor = interp_quintic3(P_scaled - P_floor);
    
    GRADIENT_COORD(0, 0, 0)
	GRADIENT_COORD(0, 0, 1)
	GRADIENT_COORD(0, 1, 0)
	GRADIENT_COORD(0, 1, 1)
	GRADIENT_COORD(1, 0, 0)
	GRADIENT_COORD(1, 0, 1)
	GRADIENT_COORD(1, 1, 0)
	GRADIENT_COORD(1, 1, 1)

    float3 P0y = mix(mix(P000, P100, P_floor.x), mix(P010, P110, P_floor.x), P_floor.y);
    float3 P1y = mix(mix(P001, P101, P_floor.x), mix(P011, P111, P_floor.x), P_floor.y);

    return P + mix(P0y, P1y, P_floor.z) * p_amplitude;
}

static inline float3 perturb_gradient_fractal(float3 P, int seed,  const int p_octaves, const float p_gain, const float p_lacunarity, const float p_amplitude, const float p_frequency)
{
    int seed_fractal = seed - 1;
    float3 P_fractal = P;
	float frequency_fractal = p_frequency;
	float amplitude_fractal = p_amplitude;
	
	P_fractal = perturb_gradient_single(P_fractal, seed_fractal, amplitude_fractal, frequency_fractal);
	
	int octave = 0;
	while (octave < p_octaves)
	{
	    frequency_fractal *= p_lacunarity;
		seed_fractal -= 1;
		amplitude_fractal *= p_gain;
		
		P_fractal = perturb_gradient_single(P_fractal, seed_fractal, amplitude_fractal, frequency_fractal);
        octave += 1;
	}

    return P_fractal;
}

static inline float3 perturb_normalize(float3 P, const float normalize_length)
{
    return P * (normalize_length * (1.0f / sqrt(dot(P, P))));
}

static inline float3 perturb_gradient_single_normalize(float3 P, int seed,  const float p_amplitude, const float p_frequency, const float normalize_length)
{
    float3 P_gradient = perturb_gradient_single(P, seed, p_amplitude, p_frequency);
    return perturb_normalize(P_gradient, normalize_length);
}

static inline float3 perturb_gradient_fractal_normalize(float3 P, int seed,  const int p_octaves, const float p_gain, const float p_lacunarity, const float p_amplitude, const float p_frequency, const float normalize_length)
{
    float3 P_gradient = perturb_gradient_fractal(P, seed, p_octaves, p_gain, p_lacunarity, p_amplitude, p_frequency);
    return perturb_normalize(P_gradient, normalize_length);
}


/* Macro definition for different fractal types.
All macros expect the following parameters to be defined:
- P : float3 = the position float3 
- seed : int = the noise_func seed
- octaves : int = fractal octaves
- gain : float = gain for each octave (also known as roughness)
- lacunarity : float = the lacunarity of each fractal octave
*/

#define FBM_SINGLE(NOISE_FUNC)\
    int seed_fractal = seed;\
    float3 P_fractal = P;\
    float amp_fractal = 1.0f;\
    int octave = 0;\
    float result = NOISE_FUNC(P_fractal, seed_fractal);\
    \
    while(octave++ < octaves)\
    {\
        P_fractal *= lacunarity;\
        seed_fractal += 1;\
        amp_fractal *= gain;\
        result += (float)(NOISE_FUNC(P_fractal, seed_fractal)) * amp_fractal;\
    }\
    result *= fractal_bound;

#define BILLOW_SINGLE(NOISE_FUNC)\
	int seed_fractal = seed;\
    float3 P_fractal = P;\
	float result = fabs((float)(NOISE_FUNC(P_fractal, seed_fractal))) * 2.0f - 1.0f;\
	float amp_fractal = 1.0f;\
	int octave = 0;\
	\
	while (octave++ < octaves)\
	{\
		P_fractal *= lacunarity;\
		seed_fractal += 1;\
		amp_fractal *= gain;\
		result += (fabs((float)(NOISE_FUNC(P_fractal, seed_fractal))) * 2.0f - 1.0f) * amp_fractal;\
    }\
	result *= fractal_bound;

#define RIGIDMULTI_SINGLE(NOISE_FUNC)\
    int seed_fractal = seed;\
    float3 P_fractal = P;\
    float amp_fractal = 1.0f;\
    int octave = 0;\
    float result = 1.0f - fabs((float)(NOISE_FUNC(P_fractal, seed_fractal)));\
	\
	while (octave++ < octaves)\
	{\
		P_fractal *= lacunarity;\
		seed_fractal += 1;\
		amp_fractal *= gain;\
		result += (1.0f - fabs((float)(NOISE_FUNC(P_fractal, seed_fractal)))) * amp_fractal;\
    }

/* All fractal noise definition for the provided noise function.*/

#define FRACTAL_NOISE_SINGLE(NOISE_FUNC)\
    static inline float fbm_##NOISE_FUNC(float3 P, int seed, const int octaves, const float lacunarity, const float gain, const float fractal_bound)\
    {\
        FBM_SINGLE(NOISE_FUNC)\
        return result;\
    }\
    \
    static inline float billow_##NOISE_FUNC(float3 P, int seed, const int octaves, const float lacunarity, const float gain, const float fractal_bound)\
    {\
        BILLOW_SINGLE(NOISE_FUNC)\
        return result;\
    }\
    \
    static inline float rigidmulti_##NOISE_FUNC(float3 P, int seed, const int octaves, const float lacunarity, const float gain, const float fractal_bound)\
    {\
        RIGIDMULTI_SINGLE(NOISE_FUNC)\
        return result;\
    }

#define euclidean_DISTANCE(position) dot(position, position)
#define manhattan_DISTANCE(position) sum3(fabs(position))
#define natural_DISTANCE(position) (euclidean_DISTANCE(position) + manhattan_DISTANCE(position))

#define distance2_RETURN(_distance, _distance2) (_distance2)
#define distance2_add_RETURN(_distance, _distance2) (_distance + _distance2)
#define distance2_sub_RETURN(_distance, _distance2) (_distance2 - _distance)
#define distance2_mul_RETURN(_distance, _distance2) (_distance * _distance2)
#define distance2_div_RETURN(_distance, _distance2) (_distance / _distance2)

#define CELLULAR_VALUE_SINGLE(DISTANCE_FUNC) \
    static inline float cellular_value_##DISTANCE_FUNC(float3 P, int seed,  const float jitter) \
    {\
        float distance_closest = 999999.0f; \
        float final_value = 0.0f; \
        \
        float3 P_base_coords_i = rint(P) - 1; \
        float3 P_base_coords_f = P_base_coords_i - P; \
        P_base_coords_i *= PRIME_VECTOR; \
        \
        float3 P_coords_i = P_base_coords_i; \
        float3 P_coords_f = P_base_coords_f; \
        for (int xi=0; xi<3; xi++) \
        { \
            P_coords_f.y = P_base_coords_f.y; \
            P_coords_i.y = P_base_coords_i.y; \
            for (int yi=0; yi<3; yi++) \
            { \
                P_coords_f.z = P_base_coords_f.z; \
                P_coords_i.z = P_base_coords_i.z; \
                for (int zi=0; zi<3; zi++) \
                { \
                    int hash = hash_vector_hb(P_coords_i, seed); \
                    float3 P_distance = (float3)( \
                        (float)(hash & BIT10_MASK), \
                        (float)((hash >> 10) & BIT10_MASK), \
                        (float)((hash >> 20) & BIT10_MASK) \
                    ) - N_511_5; \
                    \
                    P_distance = P_coords_f + (P_distance * jitter * (1.0f / sqrt(dot(P_distance, P_distance)))); \
                    \
                    float cell_value = (float)(hash) * HASH_TO_FLOAT; \
                    float distance_cell = DISTANCE_FUNC##_DISTANCE(P_distance); \
                    \
                    if(distance_cell < distance_closest) \
                    { \
                        final_value = cell_value; \
                        distance_closest = distance_cell; \
                    } \
                    \
                    P_coords_f.z += 1.0f; \
                    P_coords_i.z += Z_PRIME; \
                } \
                P_coords_f.y += 1.0f; \
                P_coords_i.y += Y_PRIME; \
            } \
            P_coords_f.x += 1.0f; \
            P_coords_i.x += X_PRIME; \
        } \
        \
        return final_value; \
    }

#define CELLULAR_DISTANCE_SINGLE(DISTANCE_FUNC) \
    static inline float cellular_distance_##DISTANCE_FUNC(float3 P, int seed, const float jitter) \
    {\
        float distance_closest = 999999.0f; \
        float final_value = 0.0f; \
        \
        float3 P_base_coords_i = rint(P) - 1; \
        float3 P_base_coords_f = P_base_coords_i - P; \
        P_base_coords_i *= PRIME_VECTOR; \
        \
        float3 P_coords_i = P_base_coords_i; \
        float3 P_coords_f = P_base_coords_f; \
        for (int xi=0; xi<3; xi++) \
        { \
            P_coords_f.y = P_base_coords_f.y; \
            P_coords_i.y = P_base_coords_i.y; \
            for (int yi=0; yi<3; yi++) \
            { \
                P_coords_f.z = P_base_coords_f.z; \
                P_coords_i.z = P_base_coords_i.z; \
                for (int zi=0; zi<3; zi++) \
                { \
                    int hash = hash_vector_hb(P_coords_i, seed); \
                    float3 P_distance = (float3)( \
                        (float)(hash & BIT10_MASK), \
                        (float)((hash >> 10) & BIT10_MASK), \
                        (float)((hash >> 20) & BIT10_MASK) \
                    ) - N_511_5; \
                    \
                    P_distance = P_coords_f + (P_distance * jitter * (1.0f / sqrt(dot(P_distance, P_distance)))); \
                    \
                    distance_closest = min(distance_closest, DISTANCE_FUNC##_DISTANCE(P_distance)); \
                    \
                    P_coords_f.z += 1.0f; \
                    P_coords_i.z += Z_PRIME; \
                } \
                P_coords_f.y += 1.0f; \
                P_coords_i.y += Y_PRIME; \
            } \
            P_coords_f.x += 1.0f; \
            P_coords_i.x += X_PRIME; \
        } \
        \
        return distance_closest; \
    }

#define CELLULAR_DISTANCE2(DISTANCE_FUNC, RETURN_FUNC) \
    static inline float cellular_distance_##DISTANCE_FUNC##_##RETURN_FUNC(float3 P, int seed, const float jitter, const float2 distance_index) \
    {\
        float distance_closest[] = {999999.0f, 999999.0f, 999999.0f, 999999.0f}; \
        float final_value = 0.0f; \
        \
        float3 P_base_coords_i = rint(P) - 1; \
        float3 P_base_coords_f = P_base_coords_i - P; \
        P_base_coords_i *= PRIME_VECTOR; \
        \
        float3 P_coords_i = P_base_coords_i; \
        float3 P_coords_f = P_base_coords_f; \
        for (int xi=0; xi<3; xi++) \
        { \
            P_coords_f.y = P_base_coords_f.y; \
            P_coords_i.y = P_base_coords_i.y; \
            for (int yi=0; yi<3; yi++) \
            { \
                P_coords_f.z = P_base_coords_f.z; \
                P_coords_i.z = P_base_coords_i.z; \
                for (int zi=0; zi<3; zi++) \
                { \
                    int hash = hash_vector_hb(P_coords_i, seed); \
                    float3 P_distance = (float3)( \
                        (float)(hash & BIT10_MASK), \
                        (float)((hash >> 10) & BIT10_MASK), \
                        (float)((hash >> 20) & BIT10_MASK) \
                    ) - N_511_5; \
                    \
                    P_distance = P_coords_f + (P_distance * jitter * (1.0f / sqrt(dot(P_distance, P_distance)))); \
                    \
                    float cell_distance = DISTANCE_FUNC##_DISTANCE(P_distance); \
                    for(int i=(int)(distance_index.y); i>0; i--) { distance_closest[i] = max(min(distance_closest[i], cell_distance), distance_closest[i-1]); } \
                    \
                    distance_closest[0] = min(distance_closest[0], cell_distance); \
                    \
                    P_coords_f.z += 1.0f; \
                    P_coords_i.z += Z_PRIME; \
                } \
                P_coords_f.y += 1.0f; \
                P_coords_i.y += Y_PRIME; \
            } \
            P_coords_f.x += 1.0f; \
            P_coords_i.x += X_PRIME; \
        } \
        \
        return RETURN_FUNC##_RETURN(distance_closest[(int)(distance_index.x)], distance_closest[(int)(distance_index.y)]); \
    }

#define CELLULAR_DISTANCE2CAVE_SINGLE(DISTANCE_FUNC) \
    static inline float cellular_distance_cave_##DISTANCE_FUNC(float3 P, int seed, const float jitter, const float2 distance_index) \
    { \
        float3 P_local = P; \
        int seed_local = seed; \
        float c0 = cellular_distance_##DISTANCE_FUNC##_distance2_div(P_local, seed_local, jitter, distance_index); \
        \
        P_local += 0.5f; \
        seed_local += 1; \
        \
        float c1 = cellular_distance_##DISTANCE_FUNC##_distance2_div(P_local, seed_local, jitter, distance_index); \
        \
        return min(c0, c1);\
    }


#define CELLULAR_LOOKUP(DISTANCE_FUNC, NOISE_FUNC) \
    static inline float cellular_value_##NOISE_FUNC##_##DISTANCE_FUNC(float3 P, int seed, const float jitter, const float frequency)\
    {\
        float distance_closest = 999999.0f; \
        float3 P_cell = P; \
        \
        float3 P_base_coords_i = rint(P) - 1; \
        float3 P_base_coords_f = P_base_coords_i; \
        P_base_coords_i *= PRIME_VECTOR; \
        \
        float3 P_coords_i = P_base_coords_i; \
        float3 P_coords_f = P_base_coords_f; \
        float3 P_local = P_coords_f - P; \
        for (int xi=0; xi<3; xi++) \
        { \
            P_local.x = P_coords_f.x - P.x; \
            P_coords_f.y = P_base_coords_f.y; \
            P_coords_i.y = P_base_coords_i.y; \
            for (int yi=0; yi<3; yi++) \
            { \
                P_local.y = P_coords_f.y - P.y; \
                P_coords_f.z = P_base_coords_f.z; \
                P_coords_i.z = P_base_coords_i.z; \
                for (int zi=0; zi<3; zi++) \
                { \
                    P_local.z = P_coords_f.z - P.z; \
                    int hash = hash_vector_hb(P_coords_i, seed); \
                    float3 P_distance = (float3)( \
                        (float)(hash & BIT10_MASK), \
                        (float)((hash >> 10) & BIT10_MASK), \
                        (float)((hash >> 20) & BIT10_MASK) \
                    ) - N_511_5; \
                    \
                    float inverse_magnitude = jitter * (1.0f / sqrt(dot(P_distance, P_distance))); \
                    float3 P_cell_new = P_distance * inverse_magnitude; \
                    P_distance = P_cell_new + P_local; \
                    P_cell_new += P_coords_f; \
                    \
                    float distance_cell = DISTANCE_FUNC##_DISTANCE(P_distance); \
                    \
                    if(distance_cell < distance_closest) \
                    { \
                        distance_closest = distance_cell; \
                        P_cell = P_cell_new; \
                    } \
                    \
                    P_coords_f.z += 1.0f; \
                    P_coords_i.z += Z_PRIME; \
                } \
                P_coords_f.y += 1.0f; \
                P_coords_i.y += Y_PRIME; \
            } \
            P_coords_f.x += 1.0f; \
            P_coords_i.x += X_PRIME; \
        } \
        \
        P_cell *= frequency; \
        \
        return NOISE_FUNC(P_cell, seed); \
    }

#define CELLULAR_FRACTAL_LOOKUP(DISTANCE_FUNC, NOISE_FUNC, FRACTAL_TYPE) \
    static inline float cellular_value_##NOISE_FUNC##_##DISTANCE_FUNC##_##FRACTAL_TYPE(float3 P, int seed, const float jitter, const float frequency, const int octaves, const float lacunarity, const float gain, const float fractal_bound)\
    {\
        float distance_closest = 999999.0f; \
        float3 P_cell; \
        \
        float3 P_base_coords_i = rint(P) - 1; \
        float3 P_base_coords_f = P_base_coords_i; \
        P_base_coords_i *= PRIME_VECTOR; \
        \
        float3 P_coords_i = P_base_coords_i; \
        float3 P_coords_f = P_base_coords_f; \
        float3 P_local = P_coords_f - P; \
        for (int xi=0; xi<3; xi++) \
        { \
            P_local.x = P_coords_f.x - P.x; \
            P_coords_f.y = P_base_coords_f.y; \
            P_coords_i.y = P_base_coords_i.y; \
            for (int yi=0; yi<3; yi++) \
            { \
                P_local.y = P_coords_f.y - P.y; \
                P_coords_f.z = P_base_coords_f.z; \
                P_coords_i.z = P_base_coords_i.z; \
                for (int zi=0; zi<3; zi++) \
                { \
                    P_local.z = P_coords_f.z - P.z; \
                    int hash = hash_vector_hb(P_coords_i, seed); \
                    float3 P_distance = (float3)( \
                        (float)(hash & BIT10_MASK), \
                        (float)((hash >> 10) & BIT10_MASK), \
                        (float)((hash >> 20) & BIT10_MASK) \
                    ) - N_511_5; \
                    \
                    float inverse_magnitude = jitter * (1.0f / sqrt(dot(P_distance, P_distance))); \
                    float3 P_cell_new = P_distance * inverse_magnitude; \
                    P_distance = P_cell_new + P_local; \
                    P_cell_new += P_coords_f; \
                    \
                    float distance_cell = DISTANCE_FUNC##_DISTANCE(P_distance); \
                    \
                    if(distance_cell < distance_closest) \
                    { \
                        distance_closest = distance_cell; \
                        P_cell = P_cell_new; \
                    } \
                    \
                    P_coords_f.z += 1.0f; \
                    P_coords_i.z += Z_PRIME; \
                } \
                P_coords_f.y += 1.0f; \
                P_coords_i.y += Y_PRIME; \
            } \
            P_coords_f.x += 1.0f; \
            P_coords_i.x += X_PRIME; \
        } \
        \
        P_cell *= frequency; \
        \
        return FRACTAL_TYPE##_##NOISE_FUNC(P_cell, seed, octaves, lacunarity, gain, fractal_bound); \
    }

#define CELLULAR_DISTANCE2_SINGLE(RETURN_FUNC) \
    CELLULAR_DISTANCE2(euclidean, RETURN_FUNC) \
    CELLULAR_DISTANCE2(manhattan, RETURN_FUNC) \
    CELLULAR_DISTANCE2(natural, RETURN_FUNC)

#define CELLULAR_LOOKUP_SINGLE(NOISE_FUNC) \
    CELLULAR_LOOKUP(euclidean, NOISE_FUNC) \
    CELLULAR_LOOKUP(manhattan, NOISE_FUNC) \
    CELLULAR_LOOKUP(natural, NOISE_FUNC)

#define CELLULAR_FRACTAL_LOOKUP_DISTANCE(NOISE_FUNC, FRACTAL_TYPE) \
    CELLULAR_FRACTAL_LOOKUP(euclidean, NOISE_FUNC, FRACTAL_TYPE) \
    CELLULAR_FRACTAL_LOOKUP(manhattan, NOISE_FUNC, FRACTAL_TYPE) \
    CELLULAR_FRACTAL_LOOKUP(natural, NOISE_FUNC, FRACTAL_TYPE)

#define CELLULAR_FRACTAL_LOOKUP_SINGLE(NOISE_FUNC) \
    CELLULAR_FRACTAL_LOOKUP_DISTANCE(NOISE_FUNC, fbm) \
    CELLULAR_FRACTAL_LOOKUP_DISTANCE(NOISE_FUNC, billow) \
    CELLULAR_FRACTAL_LOOKUP_DISTANCE(NOISE_FUNC, rigidmulti)

/* API wrappers for noise_func calls, 1D and 3D version */

static inline float value(float3 P, int seed)
{
    return value_single(P, seed);
}

static inline float3 value3(float3 P, int seed)
{
    float v = value_single(P, seed);
    return (float3)(v, v, v);
}

static inline float perlin(float3 P, int seed)
{
    return perlin_single(P, seed);
}

static inline float3 perlin3(float3 P, int seed)
{
    float v = perlin_single(P, seed);
    return (float3)(v, v, v);
}

static inline float simplex(float3 P, int seed)
{
    return simplex_single(P, seed);
}

static inline float3 simplex3(float3 P, int seed)
{
    float v = simplex_single(P, seed);
    return (float3)(v, v, v);
}

static inline float white_noise(float3 P, int seed)
{
    return white_noise_single(P, seed);
}

static inline float3 white_noise3(float3 P, int seed)
{
    float v = white_noise_single(P, seed);
    return (float3)(v, v, v);
}

static inline float cubic(float3 P, int seed)
{
    return cubic_single(P, seed);
}

static inline float3 cubic3(float3 P, int seed)
{
    float v = cubic_single(P, seed);
    return (float3)(v, v, v);
}

FRACTAL_NOISE_SINGLE(value)
FRACTAL_NOISE_SINGLE(perlin)
FRACTAL_NOISE_SINGLE(simplex)
FRACTAL_NOISE_SINGLE(white_noise)
FRACTAL_NOISE_SINGLE(cubic)

CELLULAR_VALUE_SINGLE(euclidean)
CELLULAR_VALUE_SINGLE(manhattan)
CELLULAR_VALUE_SINGLE(natural)

CELLULAR_DISTANCE_SINGLE(euclidean)
CELLULAR_DISTANCE_SINGLE(manhattan)
CELLULAR_DISTANCE_SINGLE(natural)

CELLULAR_DISTANCE2_SINGLE(distance2)
CELLULAR_DISTANCE2_SINGLE(distance2_add)
CELLULAR_DISTANCE2_SINGLE(distance2_sub)
CELLULAR_DISTANCE2_SINGLE(distance2_mul)
CELLULAR_DISTANCE2_SINGLE(distance2_div)

CELLULAR_DISTANCE2CAVE_SINGLE(euclidean)
CELLULAR_DISTANCE2CAVE_SINGLE(manhattan)
CELLULAR_DISTANCE2CAVE_SINGLE(natural)

CELLULAR_LOOKUP_SINGLE(value)
CELLULAR_LOOKUP_SINGLE(perlin)
CELLULAR_LOOKUP_SINGLE(simplex)
CELLULAR_LOOKUP_SINGLE(white_noise)
CELLULAR_LOOKUP_SINGLE(cubic) 

CELLULAR_FRACTAL_LOOKUP_SINGLE(value)
CELLULAR_FRACTAL_LOOKUP_SINGLE(perlin)
CELLULAR_FRACTAL_LOOKUP_SINGLE(simplex)
CELLULAR_FRACTAL_LOOKUP_SINGLE(white_noise)
CELLULAR_FRACTAL_LOOKUP_SINGLE(cubic)

#endif
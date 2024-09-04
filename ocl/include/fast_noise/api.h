#ifndef __fast_noise_api__
#define __fast_noise_api__

#include "fast_noise/internal.h"

#define NTYPE_VALUE 0
#define NTYPE_VALUE_FRACTAL 1
#define NTYPE_PERLIN 2
#define NTYPE_PERLIN_FRACTAL 3  
#define NTYPE_SIMPLEX 4
#define NTYPE_SIMPLEX_FRACTAL 5
#define NTYPE_WHITE_NOISE 6
#define NTYPE_CELLULAR 7
#define NTYPE_CUBIC 8
#define NTYPE_CUBIC_FRACTAL 9


#define FTYPE_FBM 0
#define FTYPE_BILLOW 1
#define FTYPE_RIGID_MULTI 2


#define DTYPE_EUCLIDEAN 0
#define DTYPE_MANHATTAN 1
#define DTYPE_NATURAL 2


#define CTYPE_CELL_VALUE 0
#define CTYPE_DISTANCE 1
#define CTYPE_DISTANCE_2 2
#define CTYPE_DISTANCE_2_ADD 3
#define CTYPE_DISTANCE_2_SUB 4
#define CTYPE_DISTANCE_2_MUL 5
#define CTYPE_DISTANCE_2_DIV 6
#define CTYPE_NOISE_LOOKUP 7
#define CTYPE_DISTANCE_2_CAVE 8

  
#define PTYPE_NONE 0
#define PTYPE_GRADIENT 1
#define PTYPE_GRADIENT_FRACTAL 2
#define PTYPE_NORMALIZE 3
#define PTYPE_GRADIENT_NORMALIZE 4
#define PTYPE_GRADIENT_FRACTAL_NORMALIZE 5

typedef struct {
    int seed;
    float frequency;
    int noise_type;
    float3 axis_scale;
    int octaves;
    float lacunarity;
    float gain;
    int fractal_type;
    int cellular_return_type;
    int cellular_distance_metric;
    int cellular_lookup_type;
    float cellular_frequency;
    float2 cellular_distance_index;
    float cellular_jitter;
    int perturbation_type;
    float perturbation_frequency;
    float perturbation_amplitude;
    int perturbation_octaves;
    float perturbation_lacunarity;
    float perturbation_gain;                     
    float perturbation_normalized_length; 
} Noise; 

static inline void noise_init(Noise * settings)
{
    settings->seed = 1337;
    settings->frequency = 0.001953125f; // -> 1.0f / 512.0f
    settings->noise_type = NTYPE_SIMPLEX;
    settings->axis_scale = (float3)(1.0f);
    settings->octaves = 1;
    settings->lacunarity = 2.01234f;
    settings->gain = 1.0f;
    settings->fractal_type = FTYPE_FBM;
    settings->cellular_return_type = CTYPE_DISTANCE;
    settings->cellular_distance_metric = DTYPE_EUCLIDEAN;
    settings->cellular_lookup_type = NTYPE_PERLIN;
    settings->cellular_frequency = 0.2f;
    settings->cellular_distance_index = (float2)(0.0f, 1.0f);
    settings->cellular_jitter = 0.45f;
    settings->perturbation_type = PTYPE_NONE;
    settings->perturbation_frequency = 0.45f;
    settings->perturbation_amplitude = 1.0f;
    settings->perturbation_octaves = 1;
    settings->perturbation_lacunarity = 2.01234f;
    settings->perturbation_gain = 1.0f;                     
    settings->perturbation_normalized_length = 1000000.0f;   
}

static inline float compute_fractal_bound(const int octaves, const float gain)
{
    float amplitude = gain;
    float f_amplitude = 1.0f;
    for (int i=1; i<octaves; i++)
    {
        f_amplitude += amplitude;
        amplitude *= gain;
    }

    return 1.0f / f_amplitude;
}

static inline float get_noise(Noise * settings, const float3 P)
{
    float n = 0.0f;
    float3 P_noise = P;

    float fractal_bound = compute_fractal_bound(settings->octaves, settings->gain);

    // Apply gradient perturbation if needed
    if(settings->perturbation_type != PTYPE_NONE)
    { 
        float normalized_length = settings->perturbation_normalized_length * settings->frequency;
        settings->perturbation_amplitude *= fractal_bound;
        
                if(settings->perturbation_type == PTYPE_GRADIENT)
        {
            P_noise = perturb_gradient_single(P_noise, settings->seed - 1, settings->perturbation_amplitude, settings->perturbation_frequency);
        }
        else if(settings->perturbation_type == PTYPE_GRADIENT_FRACTAL)
        {
            P_noise = perturb_gradient_fractal(P_noise, settings->seed - 1, settings->perturbation_octaves, settings->perturbation_gain, settings->perturbation_lacunarity, settings->perturbation_amplitude, settings->perturbation_frequency);
        }
        else if(settings->perturbation_type == PTYPE_NORMALIZE)
        {
            P_noise = perturb_normalize(P_noise, normalized_length);
        }
        else if(settings->perturbation_type == PTYPE_GRADIENT_NORMALIZE)
        {
            P_noise = perturb_gradient_single_normalize(P_noise, settings->seed - 1, settings->perturbation_amplitude, settings->perturbation_frequency, normalized_length);
        }
        else if(settings->perturbation_type == PTYPE_GRADIENT_FRACTAL_NORMALIZE)
        {
            P_noise = perturb_gradient_fractal_normalize(P_noise, settings->seed - 1, settings->perturbation_octaves, settings->perturbation_gain, settings->perturbation_lacunarity, settings->perturbation_amplitude, settings->perturbation_frequency, normalized_length);
        } 
    }  

    P_noise *= settings->axis_scale * settings->frequency;

    // Apply desired noise
            if(settings->noise_type == NTYPE_VALUE) 
    {
        n = value(P_noise, settings->seed);
    }
    else if(settings->noise_type == NTYPE_VALUE_FRACTAL)
    {
                if(settings->fractal_type == FTYPE_FBM) { n = fbm_value(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_BILLOW) { n = billow_value(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_value(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
    }
    else if(settings->noise_type == NTYPE_PERLIN)
    {
        n = perlin(P_noise, settings->seed);
    }  
    else if(settings->noise_type == NTYPE_PERLIN_FRACTAL)
    {
                if(settings->fractal_type == FTYPE_FBM) { n = fbm_perlin(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_BILLOW) { n = billow_perlin(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_perlin(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
    } 
    else if(settings->noise_type == NTYPE_SIMPLEX)
    {
        n = simplex(P_noise, settings->seed);
    } 
    else if(settings->noise_type == NTYPE_SIMPLEX_FRACTAL)
    {
                if(settings->fractal_type == FTYPE_FBM) { n = fbm_simplex(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_BILLOW) { n = billow_simplex(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_simplex(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
    }
    else if(settings->noise_type == NTYPE_WHITE_NOISE)
    {
        n = white_noise(P_noise, settings->seed);
    }  
    else if(settings->noise_type == NTYPE_CELLULAR)
    {
                if(settings->cellular_return_type == CTYPE_CELL_VALUE)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_euclidean(P_noise, settings->seed, settings->cellular_jitter); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_manhattan(P_noise, settings->seed, settings->cellular_jitter); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_natural  (P_noise, settings->seed, settings->cellular_jitter); }
            
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean(P_noise, settings->seed, settings->cellular_jitter); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan(P_noise, settings->seed, settings->cellular_jitter); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural  (P_noise, settings->seed, settings->cellular_jitter); }
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE_2)
        { 
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE_2_ADD)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_add(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_add(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_add  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE_2_SUB)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_sub(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_sub(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_sub  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE_2_MUL)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_mul(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_mul(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_mul  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE_2_DIV)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_div(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_div(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_div  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
        }
        else if(settings->cellular_return_type == CTYPE_NOISE_LOOKUP)
        {   
                    if(settings->cellular_lookup_type == NTYPE_VALUE)  
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_value_euclidean(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_value_manhattan(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_value_natural  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
            }
            else if(settings->cellular_lookup_type == NTYPE_VALUE_FRACTAL)
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_value_euclidean_fbm        (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_value_euclidean_billow     (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_value_euclidean_rigidmulti (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_value_manhattan_fbm        (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_value_manhattan_billow     (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_value_manhattan_rigidmulti (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_value_natural_fbm          (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_value_natural_billow       (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_value_natural_rigidmulti   (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                } 
            }
            else if(settings->cellular_lookup_type == NTYPE_PERLIN) 
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_perlin_euclidean(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_perlin_manhattan(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_perlin_natural  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }  
            } 
            else if(settings->cellular_lookup_type == NTYPE_PERLIN_FRACTAL)
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_perlin_euclidean_fbm       (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_perlin_euclidean_billow    (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_perlin_euclidean_rigidmulti(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_perlin_manhattan_fbm       (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_perlin_manhattan_billow    (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_perlin_manhattan_rigidmulti(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_perlin_natural_fbm         (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_perlin_natural_billow      (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_perlin_natural_rigidmulti  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
            }
            else if(settings->cellular_lookup_type == NTYPE_SIMPLEX)
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_simplex_euclidean(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_simplex_manhattan(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_simplex_natural  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
            } 
            else if(settings->cellular_lookup_type == NTYPE_SIMPLEX_FRACTAL)
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_simplex_euclidean_fbm          (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_simplex_euclidean_billow       (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_simplex_euclidean_rigidmulti   (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_simplex_manhattan_fbm          (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_simplex_manhattan_billow       (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_simplex_manhattan_rigidmulti   (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)
                {
                            if(settings->fractal_type == FTYPE_FBM)         { n = cellular_value_simplex_natural_fbm            (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_BILLOW)      { n = cellular_value_simplex_natural_billow         (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                    else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_simplex_natural_rigidmulti     (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
                }
            }
            else if(settings->cellular_lookup_type == NTYPE_WHITE_NOISE)
            {
                        if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_white_noise_euclidean(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_white_noise_manhattan(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
                else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_white_noise_natural  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_frequency); }
            }
        }
        else if(settings->cellular_return_type == CTYPE_DISTANCE_2_CAVE)
        {
                    if(settings->cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_cave_euclidean(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_cave_manhattan(P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
            else if(settings->cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_cave_natural  (P_noise, settings->seed, settings->cellular_jitter, settings->cellular_distance_index); }
        }
    }
    else if(settings->noise_type == NTYPE_CUBIC)
    {
        n = cubic(P_noise, settings->seed);
    }
    else if(settings->noise_type == NTYPE_CUBIC_FRACTAL)
    {
                if(settings->fractal_type == FTYPE_FBM) { n = fbm_cubic(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_BILLOW) { n = billow_cubic(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
        else if(settings->fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_cubic(P_noise, settings->seed, settings->octaves, settings->lacunarity, settings->gain, fractal_bound); }
    }  

    return n;
}   

#endif
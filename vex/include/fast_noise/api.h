#ifndef __fast_noise_api__
#define __fast_noise_api__

/* This module contains the public API of the vex porting of FastNoise library.
The noise definition is handled through a 'noise_settings' struct containing all
required parameters. This struct exposes a 'get_noise' function that performs
required checks and calls the appropriate implementation, contained in the
'internal' module included below */

#include "internal.h" 

/*#define NTYPE_VALUE           (0)
#define NTYPE_VALUE_FRACTAL   (1)
#define NTYPE_PERLIN          (2)
#define NTYPE_PERLIN_FRACTAL  (3)  
#define NTYPE_SIMPLEX         (4)
#define NTYPE_SIMPLEX_FRACTAL (5)
#define NTYPE_WHITE_NOISE     (6)
#define NTYPE_CELLULAR        (7)
#define NTYPE_CUBIC           (8)
#define NTYPE_CUBIC_FRACTAL   (9)*/

/* Enum for different noise type */
#define NTYPE_DO_FRACTAL  (1)
#define NTYPE_VALUE       (2)
#define NTYPE_PERLIN      (4)
#define NTYPE_SIMPLEX     (8)
#define NTYPE_CUBIC       (16)
#define NTYPE_WHITE_NOISE (32)
#define NTYPE_CELLULAR    (64)

/* Enum for different fractal type */
#define FTYPE_FBM         (0)
#define FTYPE_BILLOW      (1)
#define FTYPE_RIGID_MULTI (2)

/* Enum for different distance function in cellular mode */
#define DTYPE_EUCLIDEAN (0)
#define DTYPE_MANHATTAN (1)
#define DTYPE_NATURAL   (2)

/* Enum for different distance combination in cellular mode */
#define CTYPE_CELL_VALUE      (0)
#define CTYPE_NOISE_LOOKUP    (1)
#define CTYPE_DISTANCE        (2)
#define CTYPE_DISTANCE_2      (3)
#define CTYPE_DISTANCE_2_ADD  (4)
#define CTYPE_DISTANCE_2_SUB  (5)
#define CTYPE_DISTANCE_2_MUL  (6)
#define CTYPE_DISTANCE_2_DIV  (7)
#define CTYPE_DISTANCE_2_CAVE (8) 

/*#define PTYPE_NONE                       (0)
#define PTYPE_GRADIENT                   (1)
#define PTYPE_GRADIENT_FRACTAL           (2)
#define PTYPE_NORMALIZE                  (3)
#define PTYPE_GRADIENT_NORMALIZE         (4)
#define PTYPE_GRADIENT_FRACTAL_NORMALIZE (5)*/

/* Enum for different perturbation type */
#define PTYPE_NONE             (0)
#define PTYPE_DO_NORMALIZE     (1)
#define PTYPE_GRADIENT         (2)
#define PTYPE_GRADIENT_FRACTAL (4)

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

/* Noise definition struct - this contains all required parameters for generation */
struct noise_settings {
    int seed;
    float frequency;
    int noise_type;
    vector axis_scale;
    int octaves;
    float lacunarity;
    float gain;
    int fractal_type;
    int cellular_return_type;
    int cellular_distance_metric;
    int cellular_lookup_type;
    int cellular_lookup_fractal_type;
    float cellular_lookup_frequency;
    int cellular_lookup_octaves;
    float cellular_lookup_lacunarity;
    float cellular_lookup_gain;
    vector2 cellular_distance_index;
    float cellular_jitter;
    int perturbation_type;
    float perturbation_frequency;
    float perturbation_amplitude;
    int perturbation_octaves;
    float perturbation_lacunarity;
    float perturbation_gain;
    float perturbation_normalized_length;

    void init()
    {
        seed = 1337;
        frequency = 0.001953125f; // -> 1.0f / 512.0f
        noise_type = NTYPE_SIMPLEX;
        axis_scale = set(1.0f, 1.0f, 1.0f);
        octaves = 1;
        lacunarity = 2.01234f;
        gain = 1.0f;
        fractal_type = FTYPE_FBM;
        cellular_return_type = CTYPE_DISTANCE;
        cellular_distance_metric = DTYPE_EUCLIDEAN;
        cellular_lookup_type = NTYPE_PERLIN;
        cellular_lookup_fractal_type = FTYPE_FBM;
        cellular_lookup_frequency = 0.2f;
        cellular_lookup_octaves = 1;
        cellular_lookup_lacunarity = 2.01234f;
        cellular_lookup_gain = 1.0f;
        cellular_distance_index = set(0.0f, 1.0f);
        cellular_jitter = 0.45f;
        perturbation_type = PTYPE_NONE;
        perturbation_frequency = 0.45f;
        perturbation_amplitude = 1.0f;
        perturbation_octaves = 1;
        perturbation_lacunarity = 2.01234f;
        perturbation_gain = 1.0f;
        perturbation_normalized_length = 1_000_000.0f;
    }

    string to_string()
    {
        function string noise_type_to_str(const int ntype)
        {
            string values[] = { "VALUE", "PERLIN", "SIMPLEX", "CUBIC", "WHITE_NOISE", "CELLULAR" };
            int fractal = is_flag_set(ntype, NTYPE_DO_FRACTAL);
            return sprintf("%s%s", values[int(log(ntype - fractal) / log(2)) -1], (fractal > 0 ? "_FRACTAL" : ""));
        }

        function string perturbation_type_to_str(const int pert_type)
        {
            if(pert_type == 0) { return "NONE"; }
            if(pert_type == PTYPE_DO_NORMALIZE) { return "NORMALIZE"; }
            string values[] = { "GRADIENT", "GRADIENT_FRACTAL" };
            int normalize = is_flag_set(pert_type, PTYPE_DO_NORMALIZE);
            return sprintf("%s%s", values[int(log(pert_type - normalize) / log(2)) -1], (normalize > 0 ? "_NORMALIZE" : ""));
        }

        function string fractal_type_to_str(const int ftype)
        {
            string values[] = { "NONE", "FBM", "BILLOW", "RIGID_MULTI" };
            return values[ftype];
        }

        function string cellular_return_type_to_str(const int cell_type)
        {
            string values[] = { "VALUE", "LOOKUP", "DISTANCE", "DISTANCE_2", "DISTANCE_2_ADD", "DISTANCE_2_SUB", "DISTANCE_2_MUL", "DISTANCE_2_DIV", "DISTANCE_2_CAVE" };
            return values[cell_type];
        }

        function string cellular_distance_metric_to_str(const int metric_type)
        {
            string values[] = { "EUCLIDEAN", "MANHATTAN", "NATURAL" };
            return values[metric_type];
        }

        function string base_settings()
        {
            return sprintf(
                "BASE: noise_type: %s; seed: %i; frequency: %f; axis_scale = %f\n",
                this->noise_type_to_str(noise_type), seed, frequency, axis_scale);
        }
        
        function string fractal_settings()
        {
            return sprintf(
                "FRACTAL: fractal_type: %s; octaves : %i; lacunarity: %f, gain: %f\n",
                this->fractal_type_to_str(fractal_type), octaves, lacunarity, gain);
        }

        function string perturb_settings()
        {
            return sprintf(
                "PERTURBATION: type: %s; frequency: %f; amplitude: %f; octaves: %i; lacunarity: %f; gain: %f; normalize_length: %f;\n",
                this->perturbation_type_to_str(perturbation_type), perturbation_frequency, perturbation_amplitude,
                perturbation_octaves, perturbation_lacunarity, perturbation_gain, perturbation_normalized_length
            );
        }

        function string cellular_settings()
        {
            return sprintf(
                "CELLULAR: return_type: %s; distance_metric: %s; jitter: %f; distance_idx: %i; lookup_type: %s; fractal_type: %s; frequency: %f; octaves: %i; lacunarity: %f; gain: %f;\n",
                this->cellular_return_type_to_str(cellular_return_type), this->cellular_distance_metric_to_str(cellular_distance_metric), cellular_jitter, cellular_distance_index,
                this->noise_type_to_str(cellular_lookup_type), this->fractal_type_to_str(cellular_lookup_fractal_type), cellular_lookup_frequency, cellular_lookup_octaves,
                cellular_lookup_lacunarity, cellular_lookup_gain         
            );
        }

        return sprintf("%s%s%s%s",
            this->base_settings(), this->fractal_settings(),
            this->perturb_settings(), this->cellular_settings());            
    }

    float compute_fractal_bound(const int octaves; const float gain)
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

    float get_noise(const vector P)
    {
        float n = 0.0f;
        vector P_noise = P;

        float fractal_bound = this->compute_fractal_bound(octaves, gain);

        // Apply gradient perturbation if needed
        if(perturbation_type != PTYPE_NONE)
        { 
            float normalized_length = perturbation_normalized_length * frequency;
            perturbation_amplitude *= fractal_bound;
            
            if(is_flag_set(perturbation_type, PTYPE_DO_NORMALIZE))
            {
                if(perturbation_type == PTYPE_DO_NORMALIZE)
                {
                    P_noise = perturb_normalize(P_noise, normalized_length);
                }
                if(is_flag_set(perturbation_type,  PTYPE_GRADIENT))
                {
                    P_noise = perturb_gradient_single_normalize(P_noise, seed - 1, perturbation_amplitude, perturbation_frequency, normalized_length);
                }
                else if(is_flag_set(perturbation_type, PTYPE_GRADIENT_FRACTAL))
                {
                    P_noise = perturb_gradient_fractal_normalize(P_noise, seed - 1, perturbation_octaves, perturbation_gain, perturbation_lacunarity, perturbation_amplitude, perturbation_frequency, normalized_length);
                }
            }
            else
            {
                if(is_flag_set(perturbation_type,  PTYPE_GRADIENT))
                {
                    P_noise = perturb_gradient_single(P_noise, seed - 1, perturbation_amplitude, perturbation_frequency);
                }
                else if(is_flag_set(perturbation_type, PTYPE_GRADIENT_FRACTAL))
                {
                    P_noise = perturb_gradient_fractal(P_noise, seed - 1, perturbation_octaves, perturbation_gain, perturbation_lacunarity, perturbation_amplitude, perturbation_frequency);
                }
            }
        }  

        P_noise *= axis_scale * frequency;
        if(is_flag_set(noise_type, NTYPE_WHITE_NOISE)) { n = white_noise(P_noise, seed); }
        else if(is_flag_set(noise_type, NTYPE_VALUE))
        {
            if(is_flag_set(noise_type, NTYPE_DO_FRACTAL))
            {
                     if(fractal_type == FTYPE_FBM)         { n = fbm_value(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_BILLOW)      { n = billow_value(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_value(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
            }
            else{ n = value(P_noise, seed); }
        }
        else if(is_flag_set(noise_type, NTYPE_PERLIN))
        {
            if(is_flag_set(noise_type, NTYPE_DO_FRACTAL))
            {
                     if(fractal_type == FTYPE_FBM)         { n = fbm_perlin(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_BILLOW)      { n = billow_perlin(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_perlin(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
            }
            else { n = perlin(P_noise, seed); }
        }
        else if(is_flag_set(noise_type, NTYPE_SIMPLEX))
        {
            if(is_flag_set(noise_type, NTYPE_DO_FRACTAL))
            {
                     if(fractal_type == FTYPE_FBM)         { n = fbm_simplex(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_BILLOW)      { n = billow_simplex(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_simplex(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
            }
            else { n = simplex(P_noise, seed); }
        }
        else if(is_flag_set(noise_type, NTYPE_CUBIC))
        {
            if(is_flag_set(noise_type, NTYPE_DO_FRACTAL))
            {
                     if(fractal_type == FTYPE_FBM)         { n = fbm_cubic(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_BILLOW)      { n = billow_cubic(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
                else if(fractal_type == FTYPE_RIGID_MULTI) { n = rigidmulti_cubic(P_noise, seed, octaves, lacunarity, gain, fractal_bound); }
            }
            else { n = cubic(P_noise, seed); }
        }
        else if(is_flag_set(noise_type, NTYPE_CELLULAR))
        {
                 if(cellular_return_type == CTYPE_CELL_VALUE)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_euclidean(P_noise, seed, cellular_jitter); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_manhattan(P_noise, seed, cellular_jitter); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_natural  (P_noise, seed, cellular_jitter); }
                
            }
            else if(cellular_return_type == CTYPE_DISTANCE)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean(P_noise, seed, cellular_jitter); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan(P_noise, seed, cellular_jitter); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural  (P_noise, seed, cellular_jitter); }
            }
            else if(cellular_return_type == CTYPE_DISTANCE_2)
            { 
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2  (P_noise, seed, cellular_jitter, cellular_distance_index); }
            }
            else if(cellular_return_type == CTYPE_DISTANCE_2_ADD)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_add(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_add(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_add  (P_noise, seed, cellular_jitter, cellular_distance_index); }
            }
            else if(cellular_return_type == CTYPE_DISTANCE_2_SUB)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_sub(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_sub(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_sub  (P_noise, seed, cellular_jitter, cellular_distance_index); }
            }
            else if(cellular_return_type == CTYPE_DISTANCE_2_MUL)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_mul(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_mul(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_mul  (P_noise, seed, cellular_jitter, cellular_distance_index); }
            }
            else if(cellular_return_type == CTYPE_DISTANCE_2_DIV)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_euclidean_distance2_div(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_manhattan_distance2_div(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_natural_distance2_div  (P_noise, seed, cellular_jitter, cellular_distance_index); }
            }
            else if(cellular_return_type == CTYPE_DISTANCE_2_CAVE)
            {
                     if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_distance_cave_euclidean(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_distance_cave_manhattan(P_noise, seed, cellular_jitter, cellular_distance_index); }
                else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_distance_cave_natural  (P_noise, seed, cellular_jitter, cellular_distance_index); }
            }
            else if(cellular_return_type == CTYPE_NOISE_LOOKUP)
            {   
                if(is_flag_set(cellular_lookup_type, NTYPE_WHITE_NOISE))
                {
                         if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_white_noise_euclidean(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                    else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_white_noise_manhattan(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                    else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_white_noise_natural  (P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                }
                else if(is_flag_set(cellular_lookup_type, NTYPE_VALUE))
                {
                    if(is_flag_set(cellular_lookup_type, NTYPE_DO_FRACTAL))
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_value_euclidean_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_value_euclidean_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_value_euclidean_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_value_manhattan_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_value_manhattan_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_value_manhattan_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_NATURAL)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_value_natural_fbm         (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_value_natural_billow      (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_value_natural_rigidmulti  (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        } 
                    }
                    else
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_value_euclidean(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_value_manhattan(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_value_natural  (P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                    }
                }
                else if(is_flag_set(cellular_lookup_type, NTYPE_PERLIN))
                {
                    if(is_flag_set(cellular_lookup_type, NTYPE_DO_FRACTAL))
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_perlin_euclidean_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_perlin_euclidean_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_perlin_euclidean_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_perlin_manhattan_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_perlin_manhattan_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_perlin_manhattan_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_NATURAL)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_perlin_natural_fbm         (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_perlin_natural_billow      (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_perlin_natural_rigidmulti  (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                    }
                    else
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_perlin_euclidean(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_perlin_manhattan(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_perlin_natural  (P_noise, seed, cellular_jitter, cellular_lookup_frequency); }  
                    }
                }
                else if(is_flag_set(cellular_lookup_type, NTYPE_SIMPLEX))
                {
                    if(is_flag_set(cellular_lookup_type, NTYPE_DO_FRACTAL))
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_simplex_euclidean_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_simplex_euclidean_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_simplex_euclidean_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_simplex_manhattan_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_simplex_manhattan_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_simplex_manhattan_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_NATURAL)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_simplex_natural_fbm         (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_simplex_natural_billow      (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_simplex_natural_rigidmulti  (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                    }
                    else
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_simplex_euclidean(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_simplex_manhattan(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_simplex_natural  (P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                    }
                }
                else if(is_flag_set(cellular_lookup_type, NTYPE_CUBIC))
                {
                    if(is_flag_set(cellular_lookup_type, NTYPE_DO_FRACTAL))
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_cubic_euclidean_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_cubic_euclidean_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_cubic_euclidean_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_cubic_manhattan_fbm       (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_cubic_manhattan_billow    (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_cubic_manhattan_rigidmulti(P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                        else if(cellular_distance_metric == DTYPE_NATURAL)
                        {
                                 if(cellular_lookup_fractal_type == FTYPE_FBM)         { n = cellular_value_cubic_natural_fbm         (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_BILLOW)      { n = cellular_value_cubic_natural_billow      (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                            else if(cellular_lookup_fractal_type == FTYPE_RIGID_MULTI) { n = cellular_value_cubic_natural_rigidmulti  (P_noise, seed, cellular_jitter, cellular_lookup_frequency, cellular_lookup_octaves, cellular_lookup_lacunarity, cellular_lookup_gain, fractal_bound); }
                        }
                    }
                    else
                    {
                             if(cellular_distance_metric == DTYPE_EUCLIDEAN) { n = cellular_value_cubic_euclidean(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_MANHATTAN) { n = cellular_value_cubic_manhattan(P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                        else if(cellular_distance_metric == DTYPE_NATURAL)   { n = cellular_value_cubic_natural  (P_noise, seed, cellular_jitter, cellular_lookup_frequency); }
                    }
                }
            }
        }

        return n;
    }
}
  
#endif    
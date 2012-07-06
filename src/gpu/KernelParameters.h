
/*
rmetasim_gpu: A CUDA implementation of R package "rmetasim"
Copyright (C) 2012 Thomas Fussell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __KernelParameters_H__
#define __KernelParameters_H__

#include <curand_kernel.h>

namespace rmetasim_gpu {

/* Object Declarations */
struct KernelParameters
{
    /* Basic Parameters */
    int num_individuals;
    int current_generation;
    int state;
    int habitat;

    /* Transition Matrix Locations */
    int from_state;
    int to_state;

    /* Landscape Parameters */
    int max_individuals;
    int num_classes;
    int num_habitats;
    int num_stages;
    float selfing_rate;
    bool multiple_paternity;

    /* Individual Arrays */
    int *ids;
    int *maternal_ids;
    int *paternal_ids;
    short *classes;
    short *generations_born;
    short *generations_last_changed;
    short *generations_last_reproduced;
    int *num_offspring;
    int *genotypes;

    // Shuffling Arrays
    int *shuffling_dartboard;

    /* Offspring Arrays */
    int *offspring_maternal_indices;
    int *offspring_paternal_indices;
    short *offspring_classes;
    int *offspring_genotypes;

    /* Transition Matrices */
    float *S;
    float *R;
    float *M;

    /* RNG State Array */
    curandState *rand_states;

    /* Allele Information */
//    CudaAllele *alleles;

}; /* struct KernelParameters */

} /* namespace rmetasim_gpu */

#endif /* #ifndef __KernelParameters_H__ */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

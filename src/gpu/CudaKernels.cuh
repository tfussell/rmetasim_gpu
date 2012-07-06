/*
 * Copyright 2012 Thomas Fussell
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/*! \file CudaKernels.cuh
 *  \brief Declares several simulation functions which serve safe interfaces to
 *         lower-level CUDA kernel calls.
 */

#ifndef __CudaKernels_CUH__
#define __CudaKernels_CUH__

#include "gpu/KernelParameters.h"

namespace rmetasim_gpu {
namespace cuda {

// Debugging Kernels

/*! \p PrintIndividuals is used mostly for debugging. It takes an instance of 
 *  a \p KernelParameters and prints params.num_individuals individuals using
 *  the CUDA function printf from within a kernel.
 *
 *  \param params The object representing the full CUDA simulation state
 */
void PrintIndividuals(KernelParameters &params);

// Initialization Kernels

/*! \p InitializeRandStates operates over the seeds array in parralel through a 
 *  CUDA kernel, applying the value at seeds[i] to every params.rand_states[i] 
 *  for every i in the range [0, params.max_individuals).
 *
 *  \param seeds The array containing a seed for each rand state init function
 *  \param params The object representing the full CUDA simulation state
 */
void InitializeRandStates(const unsigned long *seeds, KernelParameters &params);

// Survival Kernels

/*! \p Survive operates over every living individual in the simulation state. It
 *  is assumed that the individuals have been compacted (see 
 *  Landscape_gpu::CompactIndividuals()). Also, it is assumed that
 *  params.num_individuals has been updated to hold the current number of
 *  of living individuals. Based on the survival transition matrix, params.S, a 
 *  new state will be chosen for each individual by choosing a single element
 *  from a multinomial distribution on the column corresponding to the
 *  particular individual's class. If this is changed, the last_changed value
 *  for that individual will be set to the current simulation generation.
 *
 *  \param params The object representing the full CUDA simulation state
 */
void Survive(KernelParameters &params);

/*! \p SurviveOffspring operates over every current offspring in offspring vectors. 
 *  It is assumed that the individuals have been compacted (see 
 *  Landscape_gpu::CompactIndividuals()). Also, it is assumed that
 *  params.num_individuals has been updated to hold the current number of
 *  of offspring individuals. Based on the survival transition matrix, params.S, a 
 *  new state will be chosen for each offspring by choosing a single element
 *  from a multinomial distribution on the column corresponding to the
 *  particular individual's class.
 *
 *  \param params The object representing the full CUDA simulation state
 */
void SurviveOffspring(KernelParameters &params);

// Reproduction Kernels

/*! \p CalculateRandomNumberOfOffspring operates over every living individual in
 *  the simulation state. It is assumed that the individuals have been compacted
 *  (see Landscape_gpu::CompactIndividuals()). Also, it is assumed that
 *  params.num_individuals has been updated to hold the current number of
 *  of living individuals. Based on the reproduction matrix, params.R, and the 
 *  current "to" state, params.to_state, a number of offspring will be generated
 *  for each individual by choosing a poisson random variable with lambda equal
 *  to params.R element [from=individual's class, to=params.to_state].
 *
 *  \param params The object representing the full CUDA simulation state
 */
void CalculateRandomNumberOfOffspring(KernelParameters &params);

/*! \p SexualReproduction operates over params.max_offspring elements in
 *  params.offspring_maternal_indices and params.offspring_paternal_indices,
 *  producing a new class and genotype for each element with non negative values
 *  in both arrays. These new values are stored in params.offspring_classes and
 *  params.offspring_genotypes respectively. Importantly, this step is combined
 *  with survival such that each offspring has already had survival applied to
 *  it. Only those which remain alive after survival are added to the landscape.
 *  This saves both computation and memory, but complicates the simulation order
 *  .
 *
 *  \param params The object representing the full CUDA simulation state
 */
void SexualReproduction(KernelParameters &params);

/*! \p FindMates has two similar but different functions. If multiple_paternity
 *  is true, it will find a new value for params.offspring_paternal_indices
 *  based on params.M and choosing a father of that class randomly from the 
 *  population for that element. If multiple_paternity is false, it will instead
 *  find a new value for params.offspring_paternal_indices only for elements
 *  with params.num_offspring positive in the same manner outlined in the case
 *  of multiple_paternity.
 *
 *  \param params The object representing the full CUDA simulation state
 */
void FindMates(KernelParameters &params);

// Miscellaneous Kernels

/*! \p FillShuffleVector is used in shuffling subsets of the global Landscape. 
 *  This algorithms is known as the dartboard method of shuffling. In this case,
 *  params.offspring_maternal_indices is used to hold shuffle data. Every
 *  element in the range [0, params.num_individuals) is mapped to the larger
 *  range [0, params.max_individuals). It is assumed that num_individuals is
 *  significantly less than max_individuals. Moreover, it is assumed that
 *  params.offspring_maternal_indices has been prefilled to -1 for every value.
 *
 *  \param params The object representing the full CUDA simulation state
 */
void FillShuffleVector(KernelParameters &params);

} // end cuda
} // end rmetasim_gpu

#endif // __CudaKernels_CUH__

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

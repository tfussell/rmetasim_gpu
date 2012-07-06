
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

#ifndef __Landscape_gpu_H__
#define __Landscape_gpu_H__

#include <list>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <curand_kernel.h>

#include "Landscape.h"
#include "PackedIndividual.h"
#include "gpu/KernelParameters.h"
#include "gpu/SimulationTimer.h"
#include "gpu/CudaIndividual.h"

namespace rmetasim_gpu {

// Object Declarations
class Landscape_gpu
{
public:
    // Constructors/Destructors
    Landscape_gpu();
    ~Landscape_gpu();

    // Host<=>Device Memory Transfer
    void FromLandscape(const unsigned int rand_seed, Landscape_statistics &landscape);
    void ToLandscape(Landscape_statistics &landscape);

    // Simulation Interface
    void Simulate(const int num_iterations, const bool compress, const int by_population);
    int PopSize(const int habitat = -1);
    void Extirpate();
    void ReproduceAndSurvive();
    void Survive();
    void LambdaAdjust(const int by_population = 1);
    void HabCarry(const int habitat = -1);
    void LandCarry();
    void CarryState(const int max_size, const int state);
    void Advance();

private:

    // Initialization
    void CalculateCapacities();
    void SetTransitionMatrices();
    void InitializeRandStates();
    void InitializeKernelParameters();
    void ResizeVectors();

    // Landscape_gpu<=>Landscape Individuals Transfer
    void ExtractIndividualsFromLandscape(Landscape_statistics &landscape, std::vector<CudaIndividual> &individuals);
    void InjectIndividualsIntoLandscape(const std::list<CudaIndividual> &individuals, Landscape_statistics &landscape);

    // Getters/Setters
    void GetIndividuals(std::list<CudaIndividual> &individuals);
    void SetIndividuals(const std::vector<CudaIndividual> &individuals);

    // Private Simulation Methods
    void CalculateMaleGameteMatrix();
    void CompactIndividuals();
    void SortIndividuals();
    void ShuffleIndividuals();
    int CountIndividuals();
    int CountIndividualsInHabitat(int habitat);
    int CountIndividualsInClass(int state);

    // Debug Methods
    void PrintIndividuals(const int max_to_print = 10);
    void PrintLandscape();
    void PrintSimulationStatistics();

    // Simulation Objects
    KernelParameters kernel_params_;

    // Debug Timing
    SimulationTimer timer_;

    // PRNG
    unsigned long rand_seed_;
    thrust::device_vector<curandState> rand_states_;

    // Landscape Parameters
    bool initialized_;
    int num_habitats_;
    int num_stages_;
    int num_classes_;
    int max_size_;
    int carrying_capacity_;
    bool multiple_paternity_;
    float selfing_rate_;
    short sum_loci;
    int current_generation_;
    int next_individual_id_;

    std::vector<double> habitat_extinction_rates_;
    std::vector<int> habitat_carrying_capacities_;

    // Individual Vectors
    thrust::device_vector<int> ids_;
    thrust::device_vector<int> maternal_ids_;
    thrust::device_vector<int> paternal_ids_;
    thrust::device_vector<short> classes_;
    thrust::device_vector<short> generations_born_;
    thrust::device_vector<short> generations_last_changed_;
    thrust::device_vector<short> generations_last_reproduced_;
    thrust::device_vector<int> num_offspring_;
    thrust::device_vector<int> genotypes_;

    // Index Shuffling Vector
    thrust::device_vector<int> shuffling_dartboard_;

    // Offspring Vectors
    thrust::device_vector<int> offspring_maternal_indices_;
    thrust::device_vector<int> offspring_paternal_indices_;
    thrust::device_vector<short> offspring_classes_;
    thrust::device_vector<int> offspring_genotypes_;

    // Transition Matrices
    thrust::host_vector<float> host_S_;
    thrust::host_vector<float> host_R_;
    thrust::host_vector<float> host_M_;
    thrust::device_vector<float> S_;
    thrust::device_vector<float> R_;
    thrust::device_vector<float> M_;

}; // end class Landscape_gpu

} // end rmetasim_gpu

#endif // __Landscape_gpu_H__

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

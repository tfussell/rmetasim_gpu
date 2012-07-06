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


/*! \file CudaKernels.cu
 *  \brief Defines several simulation function interfaces as well as associated
 *         lower-level CUDA kernel calls.
 */

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <curand_kernel.h>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "gpu/CudaKernels.cuh"
#include "const.h"
#include "gpu/const_gpu.h"

namespace rmetasim_gpu {
namespace cuda {

// Macros
#define CudaCall(X) CudaCall_((X), __FUNCTION__, __FILE__, __LINE__)

// Utility Functions

/*! \p CudaCall_ is a debugging function that is used both before and after 
 *  kernel calls to ensure no errors occured. This is not used directly but
 *  rather by the macro \p CudaCall which gives error location information.
 *
 *  \param error The CUDA error upon which this function is called
 *  \param function A string holding the name of the calling function.
 *  \param file The file in which this error was generated.
 *  \param line The line in the file at which this error ocurred.
 */
__host__
void CudaCall_(cudaError error, const char *function, const char *file, int line)
{
    if(error != cudaSuccess)
    {
	const char *errorStr = cudaGetErrorString(error);
	debug_printf("CUDA Error in %s (%s:%d): %s\n", function, file, line, errorStr);
    }

    return;
}

/*! \p CalculateNumBlocks is called before each kernel invocation to calculate
 *  the upper bound of needed blocks for a given number of elements and threads.
 *
 *  \param num_items Number of items to be operated upon
 *  \param num_threads Number of threads allowed per block, multiple of 32
 *  \return The result of (int)ceil((float)num_items/num_threads)
 */
__host__
int CalculateNumBlocks(const int num_items, const int num_threads)
{
    return (num_items - 1) / num_threads + 1;
}

/*! \p log_gamma calculates the natural logarithm of the gamma function on x.
 *  This function can only be called from device contexts.
 *
 *  \param xx The value to be used in the gamma function
 *  \return The result ln(gamma(x))
 */
__device__ 
float log_gamma(const float xx) 
{
    int j;
    float x,tmp,y,ser;

    const float cof[14] = { 57.1562356658629235,     -59.5979603554754912,
			    14.1360979747417471,     -0.491913816097620199,
			    .339946499848118887e-4,  .465236289270485756e-4,
			    -.983744753048795646e-4, .158088703224912494e-3,
			    -.210264441724104883e-3, .217439618115212643e-3,
			    -.164318106536763890e-3, .844182239838527433e-4,
			    -.261908384015814087e-4, .368991826595316234e-5 };
    
    y = x = xx;
    tmp = x + 5.24218750000000000;
    tmp = (x + 0.5) * log(tmp) - tmp;
    ser = 0.999999999999997092;

    for(j = 0; j < 14; j++) 
    {
	ser += cof[j] / ++y;
    }

    return tmp + log(2.5066282746310005 * ser / x);
}


/*! \p rpoisson generates a new random value from the poisson distribution of a
 *  given lambda. This function can only be called from device contexts.
 *
 *  \param randState A pointer to an initilized CUDA curandState
 *  \param lambda The lambda value of the poisson distribution
 *  \return The Poisson random variable
 */
__device__ 
int rpoisson(curandState *randState, const float lambda) 
{
    float u,u2,v,v2,p,t,lfac,lamexp,sqlam,loglam;
    int k;
    if(lambda < 5.) 
    {
	lamexp = exp(-lambda);
	k = -1;
	t = 1.;

	do 
	{
	    ++k;
	    t *= curand_uniform(randState);
	} while (t > lamexp);
    } 
    else 
    {
	sqlam = sqrt(lambda);
	loglam = log(lambda);

	for(;;) 
	{
	    u = 0.64 * curand_uniform(randState);
	    v = -0.68 + 1.28 * curand_uniform(randState);

	    if(lambda > 13.5) 
	    {
		v2 = v * v;
		if(v >= 0.) 
		{ 
		    if(v2 > 6.5 * u * (0.64 - u) * (u + 0.2)) 
		    {
			continue; 
		    }
		}
		else 
		{
		    if(v2 > 9.6 * u * (0.66 - u) * (u + 0.07)) 
		    {
			continue;
		    }
		}
	    }

	    k = int(floor(sqlam * (v / u) + lambda + 0.5));

	    if(k < 0)
	    {
		continue;
	    }

	    u2 = u * u;

	    if(lambda > 13.5) 
	    {
		if(v >= 0.)
		{
		    if(v2 < 15.2 * u2 * (0.61 - u) * (0.8 - u))
		    {
			break;
		    }
		}
		else 
		{
		    if (v2 < 6.76 * u2 * (0.62 - u) * (1.4 - u)) 
		    {
			break;
		    }
		}
	    }

	    lfac = log_gamma(k + 1.);
	    p = sqlam * exp(-lambda + k * loglam - lfac);

	    if (u2 < p)
	    { 
		break;
	    }
	}
    }

    return k;
}


// Debugging Kernels

__global__
void PrintIndividualsKernel(const KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;
	
    if(global_index < params.num_individuals)
    {
	printf("[%d] : id=%d, mid=%d, pid=%d, cl=%d, born=%d, changed=%d lastrep=%d\n", global_index, 
	    params.ids[global_index], params.maternal_ids[global_index], params.paternal_ids[global_index], 
	    params.classes[global_index], params.generations_born[global_index], 
	    params.generations_last_changed[global_index], params.generations_last_reproduced[global_index]);
    }
}

void PrintIndividuals(KernelParameters &params)
{
    if(params.num_individuals > 0)
    {
	const int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	CudaCall(cudaGetLastError());

	PrintIndividualsKernel<<< num_blocks, num_threads >>>(params);
	
	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}


// Initialization Kernels

__global__ 
void InitializeRandStatesKernel(const unsigned long *seeds, KernelParameters params)
{
    const int global_index = blockIdx.x * blockDim.x + threadIdx.x;

    if(global_index < params.max_individuals)
    {
	unsigned long seed = seeds[global_index];
	curandState *rand_state = params.rand_states + global_index;
	curand_init(seed, 0, 0, rand_state);
    }
}

void InitializeRandStates(const unsigned long *seeds, KernelParameters &params)
{
    if(params.max_individuals > 0)
    {
	const unsigned int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const unsigned int num_blocks = CalculateNumBlocks(params.max_individuals, num_threads);

	CudaCall(cudaGetLastError());

	InitializeRandStatesKernel<<< num_blocks, num_threads >>>(seeds, params);

	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());  
    }
}


// Survival Kernels

__global__ 
void SurviveKernel(KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;

    if(global_index < params.num_individuals)
    {
	const short old_class = params.classes[global_index];

	if(old_class >= 0)
	{
	    short new_class = -1;
	    short start_index = old_class * params.num_classes;

	    float cumulative_prob = 0;
	    float rand = curand_uniform(params.rand_states + global_index);

	    for(int j = 0; j < params.num_classes; j++)
	    {
		cumulative_prob += params.S[start_index + j];

		if(new_class == -1 && rand < cumulative_prob)
		{
		    new_class = j;
		}
	    }

	    params.classes[global_index] = new_class;

	    if(new_class != old_class)
	    {
		params.generations_last_changed[global_index] = params.current_generation;
	    }
	}
    }
}

void Survive(KernelParameters &params)
{	
    if(params.num_individuals > 0)
    {
	const unsigned int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const unsigned int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	SurviveKernel<<< num_blocks, num_threads >>>(params);

	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}

__global__ 
void SurviveOffspringKernel(KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;

    if(global_index < params.num_individuals)
    {
	const short old_class = params.offspring_classes[global_index];

	if(old_class >= 0)
	{
	    short new_class = -1;
	    short start_index = old_class * params.num_classes;

	    float cumulative_prob = 0;
	    float rand = curand_uniform(params.rand_states + global_index);

	    for(int j = 0; j < params.num_classes; j++)
	    {
		cumulative_prob += params.S[start_index + j];

		if(new_class == -1 && rand < cumulative_prob)
		{
		    new_class = j;
		}
	    }

	    params.offspring_classes[global_index] = new_class;
	}
    }
}

void SurviveOffspring(KernelParameters &params)
{	
    if(params.num_individuals > 0)
    {	
	assert(params.num_individuals <= OFFSPRING_ALLOC_CHUNK_SIZE);
	
	const unsigned int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const unsigned int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	SurviveOffspringKernel<<< num_blocks, num_threads >>>(params);

	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}


// Reproduction Kernels

__global__
void CalculateRandomNumberOfOffspringKernel(KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;

    if(global_index < params.num_individuals)
    {
	int from = params.classes[global_index];

	if(from >= 0)
	{
	    float mu = params.R[from * params.num_classes + params.to_state];
	    int poisson = rpoisson(params.rand_states + global_index, mu);

	    params.num_offspring[global_index] = poisson;
	}
    }
}

void CalculateRandomNumberOfOffspring(KernelParameters &params)
{
    if(params.num_individuals > 0)
    {
	const int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	CalculateRandomNumberOfOffspringKernel<<< num_blocks, num_threads >>>(params);

	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}

__global__
void SexualReproductionKernel(KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;

    if(global_index < params.num_individuals)
    {
	short new_class = -1;
	short start_index = params.to_state * params.num_classes;
	
	float rand = curand_uniform(params.rand_states + global_index);
	float cumulative_prob = 0;

	for(int j = 0; j < params.num_classes; j++)
	{
	    cumulative_prob += params.S[start_index + j];
	    
	    if(new_class == -1 && rand < cumulative_prob)
	    {
		new_class = j;
	    }
	}

	params.offspring_classes[global_index] = params.to_state;
	    
	int maternal_index = params.offspring_maternal_indices[global_index];
	int paternal_index = -1;
	
	if(params.multiple_paternity)
	{
	    paternal_index = params.offspring_paternal_indices[global_index];
	}
	else
	{
	    paternal_index = params.offspring_paternal_indices[maternal_index];
	}

	if(maternal_index != -1 && paternal_index != -1)
	{
	    rand = curand_uniform(params.rand_states + global_index);

	    if(rand < params.selfing_rate)
	    {
		paternal_index = params.offspring_maternal_indices[global_index];
	    }
	    
	    params.maternal_ids[global_index] = params.ids[maternal_index];
	    params.paternal_ids[global_index] = params.ids[paternal_index];
	    
	    int *mother_start = params.genotypes + maternal_index;
	    int *father_start = params.genotypes + paternal_index;
	    
	    int *self_start = params.offspring_genotypes + global_index;
	    
	    for(int j = 0; j < MAXLOCI; j++)
	    {
		char which_parent_allele = curand_uniform(params.rand_states + global_index);
		
		// Mother's Gamete
		self_start[j] = (!!(which_parent_allele & 1)) * mother_start[j] + 
		    (1 - !!(which_parent_allele & 1)) * mother_start[j + 1];
		
		// Father's Gamete
		self_start[j + 1] = (!!(which_parent_allele & 2)) * father_start[j] + 
		    (1 - !!(which_parent_allele & 2)) * father_start[j + 1];
	    }
	}
	else
	{
	    printf("Error in SexualReproductionKernel(): no parent found for offspring\n");
	}
    }
}

void SexualReproduction(KernelParameters &params)
{
    if(params.num_individuals > 0)
    {
	const int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	CudaCall(cudaGetLastError());

	SexualReproductionKernel<<< num_blocks, num_threads >>>(params);

	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}

__global__
void FindMatesKernel(KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;

    if(global_index < params.num_individuals)
    {
	if(params.multiple_paternity || 
	    (!params.multiple_paternity && (params.num_offspring[global_index] >= 0)))
	{
	    const short from = params.classes[global_index];

	    short target_father_class = -1;
	    float cumulative_prob = 0;
	    float rand = curand_uniform(params.rand_states + global_index);
	    
	    for(int to = 0; to < params.num_classes; to++)
	    {
		float m_value = params.M[from + to * params.num_classes];
		cumulative_prob += m_value;
		
		if(target_father_class == -1 && rand < cumulative_prob)
		{
		    target_father_class = to;
		}
	    }

	    //printf("FindMatesKernel: thread=%d, from=%d, target=%d\n",global_index,from,target_father_class);

	    if(target_father_class != -1)
	    {
		int father_index = curand(params.rand_states + global_index) * params.num_individuals - 1;

		short father_class = params.classes[father_index];
		short num_tries = 0;
		
		while(num_tries++ < MAX_MATE_SEARCH_ATTEMPTS && target_father_class != father_class)
		{
		    father_index = curand(params.rand_states + global_index) * params.num_individuals - 1;
		    father_class = params.classes[father_index];
		}
		
		if(father_class == target_father_class)
		{
		    params.offspring_paternal_indices[global_index] = father_index;
		}
		else
		{
		    // No father in target class found, what should we do?
		    //assert(father_class == target_father_class);
		} 
	    }
	}
    }
}

void FindMates(KernelParameters &params)
{
    if(params.num_individuals > 0)
    {
	const int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	//printf("cuda::FindMates(): num_threads=%d, num_blocks=%d\n",num_threads,num_blocks);

	CudaCall(cudaGetLastError());

	FindMatesKernel<<< num_blocks, num_threads >>>(params);

	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}


// Miscellaneous Kernels

// Note: This function will freeze program execution if there are not enough open positions.
__global__
void FillShuffleVectorKernel(KernelParameters params)
{
    const int global_index = blockDim.x * blockIdx.x + threadIdx.x;
    const int emptiness_indicator = -1;

    if(global_index < params.num_individuals
	&& params.classes[global_index] == params.state)
    {
	curandState *rand_state = params.rand_states + global_index;

	int new_location = curand(rand_state) % params.max_individuals;
	bool succeeded = atomicCAS(params.shuffling_dartboard + new_location, 
	    emptiness_indicator, global_index) == emptiness_indicator;

	while(!succeeded)
	{
	    new_location = curand(rand_state) % params.max_individuals;
	    succeeded = atomicCAS(params.shuffling_dartboard + new_location, 
		emptiness_indicator, global_index) == emptiness_indicator;
	}
    }
}

void FillShuffleVector(KernelParameters &params)
{
    if(params.num_individuals > 0)
    {
	const int num_threads = MAX_CUDA_THREADS_PER_BLOCK;
	const int num_blocks = CalculateNumBlocks(params.num_individuals, num_threads);

	CudaCall(cudaGetLastError());

	FillShuffleVectorKernel<<< num_blocks, num_threads >>>(params);
	
	CudaCall(cudaGetLastError());
	CudaCall(cudaDeviceSynchronize());
    }
}

} // end cuda
} // end rmetasim_gpu

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

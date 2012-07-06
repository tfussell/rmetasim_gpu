
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

#ifndef __const_gpu_H__
#define __const_gpu_H__

// 1 -> Enable debugging outputs.
// 0 -> Disable debugging outputs.
#define RMETASIM_GPU_DEBUG 1

// Restrict to cards of compute capability >= 2.0.
#define MIN_COMPUTE_CAPABILITY_MAJOR 2
#define MIN_COMPUTE_CAPABILITY_MINOR 0

// This is defined by compute capability.
// For small, fast kernels, use the maximum number of threads per block.
#define MAX_CUDA_THREADS_PER_BLOCK 1024

// Offspring are matched to fathers of a specific class.
// Then, a father in that class is chosen randomly.
// Stop looking if not found in this many tries.
// This mostly comes into play in small populations, hopefully not an issue.
#define MAX_MATE_SEARCH_ATTEMPTS 100

// We can have up to rpoisson(lambda * pop_size)? offspring per step.
// We don't want to make arrays that big, so we only do this number at a time.
#define OFFSPRING_ALLOC_CHUNK_SIZE 100000

// Individuals are transferred to the GPU in chunks so that large regions of
// memory do not need to be allocated to store them.
#define MAX_INDIVIDUAL_TRANSFER_SIZE 100000

// We need enough room to move all offspring into landscape, but we don't
// want to waste memory. carrying_capacity_ is mulitplied by this to determine
// max_size_. Ideally, this might be lambda + some amount?
#define LANDSCAPE_VECTOR_SIZE_MULTIPLIER 2

// All debug messages are done through this interface.
// The preprocessor can then remove these messages at compile-time.
// This uses a variadic macro, supported by most modern compilers, hopefully.
#if(RMETASIM_GPU_DEBUG == 1)
  #define debug_printf(fmt, ...) std::printf(fmt, ##__VA_ARGS__)
#else
  #define debug_printf(fmt, ...) 
#endif

#endif /* ifndef __const_gpu_H__ */

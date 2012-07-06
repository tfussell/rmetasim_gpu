
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

/* Local Includes */
#include "gpu/CudaIndividual.h"

namespace rmetasim_gpu {
namespace functors {

struct kill_in_range
{
    kill_in_range(short n, short x) : min(n), max(x) { }

    __host__ __device__ short operator()(const short &x)
    {
	bool kill = (x >= min && x < max);
	short new_class = kill * -1 + (1 - kill) * x;
	return new_class;
    }

    const short min;
    const short max;
};

struct in_range
{
    in_range(short n, short x) : min(n), max(x) { }

    __host__ __device__ bool operator()(const short &x)
    {
	return x >= min && x < max;
    }

    const short min;
    const short max;
};

struct equal
{
    equal(short n) : value(n) { }

    __host__ __device__ bool operator()(const short &x)
    {
	return x == value;
    }

    const short value;
};

struct tuple_class_equal
{
    tuple_class_equal(short n) : value(n) { }

    __host__ __device__ bool operator()(const CudaIndividualTuple &x)
    {
	return x.get<3>() == value;
    }

    const short value;
};

struct non_negative
{
    __host__ __device__ bool operator()(const short &x)
    {
	return x >= 0;
    }
};

struct non_negative_int
{
    __host__ __device__ bool operator()(const int &x)
    {
	return x >= 0;
    }
};

struct tuple_class_non_negative
{
    __host__ __device__ bool operator()(const CudaIndividualTuple &x)
    {
	return x.get<3>() >= 0;
    }
};

struct positive
{
    __host__ __device__ bool operator()(const short &x)
    {
	return x > 0;
    }
};

struct positive_int
{
    __host__ __device__ bool operator()(const int &x)
    {
	return x > 0;
    }
};

struct individidual_to_tuple
{
    __host__ __device__  CudaIndividualTuple operator()(const CudaIndividual in)
    {
	return thrust::make_tuple(in.id, in.mid, in.pid, in.cl, in.gen, in.changed, in.lastrep);
    }
};

struct tuple_to_individual
{
    __host__ __device__ CudaIndividual operator()(const CudaIndividualTuple in)
    {
	CudaIndividual ind;

	ind.id = in.get<0>();
	ind.mid = in.get<1>();
	ind.pid = in.get<2>();
	ind.cl = in.get<3>();
	ind.gen = in.get<4>();
	ind.changed = in.get<5>();
	ind.lastrep = in.get<6>();

	return ind;
    }
};

} /* namespace functors */
} /* namespace rmetasim_gpu */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

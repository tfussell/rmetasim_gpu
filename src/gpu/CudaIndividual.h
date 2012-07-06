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


/*! \file CudaIndividual.h
 *  \brief Defines CudaIndividual, an object holding information on one 
 *         PackedIndividual from metasim.
 */

#ifndef __CudaIndividual_H__
#define __CudaIndividual_H__

#include "PackedIndividual.h"

namespace rmetasim_gpu {

// Typedef used to shorten function calls and allow for centralized changes.
typedef thrust::tuple<int, int, int, short, short, short, short> 
CudaIndividualTuple;

/*! A \p CudaIndividual is a POD which can exist in GPU device memory or in
 *  system memory. It is used as an intermediate between PackedIndividual and
 *  the structure of arrays (SoA) holding equivalent data used in Landscape_gpu.
 *  Some of those members are not replicated in this object. These members were
 *  not found to be necessary and were not used to save memory.
 */
struct CudaIndividual
{
    int id;
    int mid;
    int pid;
    short cl;
    short gen;
    short changed;
    short lastrep;
}; // end CudaIndividual

} // end rmetasim_gpu

#endif // __CudaIndividual_H__

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

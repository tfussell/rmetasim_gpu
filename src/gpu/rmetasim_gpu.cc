
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
#include "gpu/Landscape_gpu.cuh"
#include "gpu/rmetasim_gpu.h" // include this last as R headers within cause conflicts if included first

extern "C" {
    extern void convert_R_to_metasim(SEXP Rland, Landscape_statistics &L);
    extern SEXP convert_metasim_to_R(Landscape_statistics &L);

    SEXP iterate_landscape_gpu(SEXP numit, SEXP Rland, SEXP cmpress, SEXP bypop)
    {
	Landscape_statistics Lstats_in, Lstats_out;
  
	convert_R_to_metasim(Rland,Lstats_in);

	Lstats_in.ChooseEpoch();
	Lstats_in.ConstructDemoMatrix();

	int n = INTEGER(coerceVector(numit,INTSXP))[0];
	bool compress = INTEGER(coerceVector(cmpress,INTSXP))[0];
	bool bp = INTEGER(coerceVector(bypop, INTSXP))[0];

	rmetasim_gpu::Landscape_gpu Lgpu;

	int seed = 123456789;

	Lgpu.FromLandscape(seed, Lstats_in);
	Lgpu.Simulate(n, compress, bp);
//	Lgpu.ToLandscape(Lstats_out);

	return convert_metasim_to_R(Lstats_out);
    }

} /// end of extern "C"

/*
  ;;; Local Variables: ***
  ;;; mode: C++ ***
  ;;; minor-mode: font-lock ***
  ;;; End: ***
*/

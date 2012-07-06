
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

#include <random>

#include "Landscape.h"
#include "Landscape_gpu.cuh"

void DebugConstructLandscape(int num_habitats, int num_stages, TransMat S, TransMat R, TransMat M, double *evec, int *kvec, double *ldemo, std::vector<int> &pop_sizes, Landscape &land)
{
	land.init(num_habitats, num_stages, 0, 1, 1, 1);
	land.setepochs(1);
	land.setepochprob(0, 1);
	land.setCepoch(0);
	land.setCgen(0);
	land.setdensdep(false);
	land.setextinct(0, &evec[0]);
	land.setgens(2);
	land.sethabs(num_habitats);
	land.setk(0, &kvec[0]);
	land.setldemovector(0, ldemo);
	land.setloci();
	land.setmultp(1);
	land.setndemo(1);
	land.setnextID(1);
	land.setranddemo(0);
	land.setself(0.0);
	land.setstages(num_stages);
	land.popsizeset(pop_sizes);
	land.setS(S,0);
	land.setR(R,0);
	land.setM(M,0);

	land.ChooseEpoch();
	land.ConstructDemoMatrix();
}

int main()
{
	int num_individuals = 4200;
	int num_habitats = 2;
	int num_stages = 3;
	int num_classes = num_habitats * num_stages;
	int num_loci = 3;
	int num_generations = 100;
	bool density_dependence = false;
	bool multiple_paternity = false;
	double selfing_probability = 0.0;

	double earr[] = { 0., 0. };
	int karr[] = { 15000, 15000 };

	float Sa[] = { 0.1f,0.4f,0.5f,
				   0.0f,0.0f,0.9f,
				   0.0f,0.0f,0.87f };

	float Ra[] = { 0.0f,0.0f,0.0f,
			 	   0.0f,0.0f,0.0f,
			 	   1.1f,0.0f,0.0f };

	float Ma[] = { 0.0,0.0,0.0,
				   0.0,0.0,0.0,
				   0.0,0.0,1.0 };

	TransMat S, R, M;

	S.SetSize(num_classes);
	R.SetSize(num_classes);
	M.SetSize(num_classes);

	for(int from = 0; from < num_classes; from++)
	{
		for(int to = 0; to < num_classes; to++)
		{
			int linear_index = from * num_classes + to;

			S.SetElement(from, to, Sa[linear_index]);
			R.SetElement(from, to, Ra[linear_index]);
			M.SetElement(from, to, Ma[linear_index]);
		}
	}

	double ldemo[] = { 1. };

	std::vector<int> pop_sizes(num_classes, num_individuals / num_classes);

	Landscape land;

	DebugConstructLandscape(num_habitats, num_stages, S, R, M, earr, karr, ldemo, pop_sizes, land);

	rmetasim_gpu::Landscape_gpu tland(land, 1989);
	
	tland.ToDevice();
	tland.Simulate(50, false, false);
	tland.FromDevice();

    return 0;
}

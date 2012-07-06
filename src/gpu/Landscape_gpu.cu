
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

/* External Includes */
#include <cassert>
#include <ctime>
#include <thrust/count.h>
#include <thrust/scan.h>
#include <thrust/tuple.h>
#include <thrust/binary_search.h>

/* Local Includes */
#include "gpu/Landscape_gpu.cuh"
#include "gpu/const_gpu.h"
#include "gpu/CudaKernels.cuh"
#include "gpu/ThrustFunctors.h"
#include "gpu/CudaIndividual.h"
#include "gpu/PRNG.h"

namespace rmetasim_gpu {

// Constructors/Destructors

Landscape_gpu::Landscape_gpu() :
    timer_(),
    rand_seed_(0),
    initialized_(false),
    num_habitats_(0),
    num_stages_(0),
    num_classes_(0),
    max_size_(0),
    carrying_capacity_(0),
    multiple_paternity_(false),
    selfing_rate_(0.0),
    sum_loci(0),
    current_generation_(0),
    next_individual_id_(0)
{

}

Landscape_gpu::~Landscape_gpu()
{

}

// Host<=>Device Memory Transfer

void Landscape_gpu::FromLandscape(const unsigned int rand_seed, Landscape_statistics &landscape)
{
    timer_.BeginMethod("Initialize");

    num_habitats_ = landscape.gethabs();
    num_stages_ = landscape.getstages();
    num_classes_ = num_habitats_ * num_stages_;
    multiple_paternity_ = landscape.getmultp();
    selfing_rate_ = landscape.getself();
    current_generation_ = landscape.getCgen();
    next_individual_id_ = landscape.getnextID();

    habitat_extinction_rates_.resize(num_habitats_);
    habitat_carrying_capacities_.resize(num_habitats_);

    landscape.getextinct(0, &habitat_extinction_rates_[0]);
    landscape.getk(0, &habitat_carrying_capacities_[0]);

    host_S_.resize(num_classes_ * num_classes_);
    host_R_.resize(num_classes_ * num_classes_);
    host_M_.resize(num_classes_ * num_classes_);

    S_.resize(num_classes_ * num_classes_);
    R_.resize(num_classes_ * num_classes_);
    M_.resize(num_classes_ * num_classes_);    

    for(int from = 0; from < num_classes_; from++)
    {
	for(int to = 0; to < num_classes_; to++)
	{
	    //XXX: should this be from*num_classes or to*num_classes?
	    int linear_index = from * num_classes_ + to;
	    
	    host_S_[linear_index] = landscape.getSmatElement(0, to, from);
	    host_R_[linear_index] = landscape.getRmatElement(0, to, from);
	    host_M_[linear_index] = landscape.getMmatElement(0, to, from);
	}
    }

    SetTransitionMatrices();

    std::vector<CudaIndividual> individuals;
    ExtractIndividualsFromLandscape(landscape, individuals);
	
    CalculateCapacities();
    ResizeVectors();

    if(individuals.size() > 0)
    {
	assert(individuals.size() <= max_size_);
	SetIndividuals(individuals);
    }

    InitializeKernelParameters();
    InitializeRandStates();
    initialized_ = true;

    timer_.EndMethod();
}

void Landscape_gpu::ToLandscape(Landscape_statistics &landscape)
{
    timer_.BeginMethod("Terminate");

    std::list<CudaIndividual> individuals;
    GetIndividuals(individuals);
    InjectIndividualsIntoLandscape(individuals, landscape);

    initialized_ = false;

    timer_.EndMethod();
}

// Simulation Interface

void Landscape_gpu::Simulate(const int num_iterations, const bool compress, const int by_population)
{
//    printf("\nBeginning Simulation\n");

    for(int i = 0; i < num_iterations; i++)
    {
//	printf("Step %d\n", i);

	const int pop_before = CountIndividuals();

	if(pop_before > 0)
	{
	    Extirpate();
	    ReproduceAndSurvive();
	    LambdaAdjust(by_population);
	    LandCarry();
	    HabCarry();
	}

	Advance();

//	int pop_after = CountIndividuals();
//	int pop_delta = pop_after - pop_before;

//	printf("Population Sizes: before=%d after=%d change=%d\n", pop_before, pop_after, pop_delta);
//	printf("Step %d End\n", i);
    }

    if(compress)
    {
	Survive();
    }

    LandCarry();
    HabCarry();

//    printf("\nSimulation Completed\n");

//    PrintSimulationStatistics();
}

int Landscape_gpu::PopSize(const int habitat)
{
    return habitat >= 0 ? CountIndividualsInHabitat(habitat) : CountIndividuals();
}

void Landscape_gpu::Extirpate()
{
    timer_.BeginMethod("Extirpate");

    bool changed = false;

    for(int habitat_index = 0; habitat_index < num_habitats_; habitat_index++)
    {
	double rn = rand() / (double)RAND_MAX;

	if(habitat_extinction_rates_[habitat_index] > rn)
	{
	    changed = true;

	    int min_class = habitat_index * num_stages_;
	    int max_class = min_class + num_stages_;

	    thrust::transform(classes_.begin(), classes_.end(), classes_.begin(), functors::kill_in_range(min_class, max_class));
	}
    }

    if(changed)
    {
	CompactIndividuals();
    }

    timer_.EndMethod();
}

void Landscape_gpu::ReproduceAndSurvive()
{
    timer_.BeginMethod("ReproduceAndSurvive");

    int current_pop_size = CountIndividuals();
    int num_parents = current_pop_size;

    CalculateMaleGameteMatrix();
	
    for(int to = 0; to < num_classes_; to++)
    {
	thrust::fill(num_offspring_.begin(), num_offspring_.begin() + current_pop_size, 0);

	kernel_params_.to_state = to;
	kernel_params_.num_individuals = current_pop_size;
	
	cuda::CalculateRandomNumberOfOffspring(kernel_params_);
	
	int num_offspring = thrust::count_if(num_offspring_.begin(), 
	    num_offspring_.begin() + current_pop_size, functors::positive_int());
		
	thrust::replace_if(generations_last_reproduced_.begin(), 
	    generations_last_reproduced_.end(), num_offspring_.begin(), 
	    functors::positive(), current_generation_);   

//	printf("Reproduce(): to=%d num_offspring=%d\n",to,num_offspring);

	if(num_offspring > 0)
	{
	    int num_chunks = (num_offspring - 1) / OFFSPRING_ALLOC_CHUNK_SIZE + 1;
	    int remaining_offspring = num_offspring;

//	    printf("Using %d chunks for offspring\n",num_chunks); 
    
	    for(int i = 0; i < num_chunks; i++)
	    {
		int chunk_size = std::min(OFFSPRING_ALLOC_CHUNK_SIZE, remaining_offspring);
		int chunk_start_index = i * OFFSPRING_ALLOC_CHUNK_SIZE;	

//		printf("Chunk %d: offspring from indices %d to %d\n", i, chunk_start_index, chunk_start_index + chunk_size);

		//Find Parents
		thrust::fill(offspring_maternal_indices_.begin(),
		    offspring_maternal_indices_.end(), -1);
		thrust::fill(offspring_paternal_indices_.begin(),
		    offspring_paternal_indices_.end(), -1);	    
		thrust::fill(offspring_classes_.begin(),
		    offspring_classes_.end(), -1);

		thrust::fill(offspring_classes_.begin(),
		    offspring_classes_.begin() + chunk_size, to);
	    
		thrust::inclusive_scan(num_offspring_.begin(), num_offspring_.end(), 
		    num_offspring_.begin());

		thrust::lower_bound(num_offspring_.begin(), num_offspring_.end(), 
		    thrust::counting_iterator<int>(chunk_start_index), 
		    thrust::counting_iterator<int>(chunk_start_index + chunk_size), 
		    offspring_maternal_indices_.begin());

		kernel_params_.to_state = to;
		kernel_params_.num_individuals = chunk_size;

		//cuda::FindMates(kernel_params_);
		// End Find Parents

		cuda::SurviveOffspring(kernel_params_);

		int num_living_offspring_in_chunk = thrust::count_if(
		    offspring_classes_.begin(), offspring_classes_.end(), 
		    functors::non_negative());

//		printf("%d offspring remaining in chunk after survive.\n", num_living_offspring_in_chunk);

		// Stream compaction of offspring vectors
		thrust::copy_if(
		    thrust::make_zip_iterator(
			thrust::make_tuple(
			    offspring_maternal_indices_.begin(),
			    offspring_paternal_indices_.begin(),
			    offspring_classes_.begin())),
		    thrust::make_zip_iterator(
			thrust::make_tuple(
			    offspring_maternal_indices_.end(),
			    offspring_paternal_indices_.end(),
			    offspring_classes_.end())),
		    offspring_classes_.begin(),
		    thrust::make_zip_iterator(
			thrust::make_tuple(
			    offspring_maternal_indices_.begin(),
			    offspring_paternal_indices_.begin(),
			    offspring_classes_.begin())),
		    functors::non_negative());

		//Offspring Block
		{
		    //cuda::SexualReproduction(kernel_params_);
		}
	
		// Transfer living offspring to real landscape
		if(num_living_offspring_in_chunk > 0)
		{
		    assert(current_pop_size + num_living_offspring_in_chunk <= max_size_);
		    
		    if(next_individual_id_ + num_living_offspring_in_chunk < MAXIDS)
		    {
			thrust::copy(thrust::make_counting_iterator(next_individual_id_), 
			    thrust::make_counting_iterator(next_individual_id_ + num_living_offspring_in_chunk), 
			    ids_.begin() + current_pop_size);
			next_individual_id_ += num_living_offspring_in_chunk;
		    }
		    else
		    {
			// ID wraparound occured; Assign ids in two sets.
			// First set is range [next_individual_id_,MAXIDS).
			// Second set is in range [1, ((next_individual_id_ + num_offspring) - MAXIDS)).
			int under = MAXIDS - next_individual_id_;
			int over = (next_individual_id_ + num_living_offspring_in_chunk) - MAXIDS;
			
//			printf("ID wraparound: next=%d num_offspring=%d under=%d over=%d\n",next_individual_id_,num_living_offspring_in_chunk,under,over);
		
			thrust::copy(thrust::make_counting_iterator(next_individual_id_), 
			    thrust::make_counting_iterator(MAXIDS), 
			    ids_.begin() + current_pop_size);
			thrust::copy(thrust::make_counting_iterator(1), 
			    thrust::make_counting_iterator(over), 
			    ids_.begin() + current_pop_size + under);
		
			next_individual_id_ = over;
		    }
	    
		    thrust::copy(offspring_classes_.begin(),
			offspring_classes_.begin() + num_living_offspring_in_chunk,
			classes_.begin() + current_pop_size);

		    thrust::copy(offspring_genotypes_.begin(),
			offspring_genotypes_.begin() + num_living_offspring_in_chunk,
			genotypes_.begin() + current_pop_size);
	    
		    thrust::fill(generations_born_.begin() + current_pop_size, 
			generations_born_.begin() + current_pop_size + num_living_offspring_in_chunk,
			current_generation_);
	    
		    thrust::fill(generations_last_changed_.begin() + current_pop_size, 
			generations_last_changed_.begin() + current_pop_size + num_living_offspring_in_chunk,
			current_generation_);
	    
		    //ReproduceInitializeNewbornOffspring();

		    current_pop_size += num_living_offspring_in_chunk;

		    // Remove this in production code.
		    // Ensures that reproduction is functioning correctly.
		    int check_pop_size = CountIndividuals();
		    assert(current_pop_size == check_pop_size);
		}
	    }
	}
    }

    //Survive Parents Only
    {
	kernel_params_.num_individuals = num_parents;
	kernel_params_.current_generation = current_generation_;
	
	cuda::Survive(kernel_params_);
	
	CompactIndividuals();
    }

    timer_.EndMethod();
}

void Landscape_gpu::Survive()
{
    timer_.BeginMethod("Survive");

    kernel_params_.current_generation = current_generation_;

    cuda::Survive(kernel_params_);
    CompactIndividuals();

    timer_.EndMethod();
}


void Landscape_gpu::LambdaAdjust(const int by_population)
{
    timer_.BeginMethod("LambdaAdjust");

    int i, j, k, l, bigto, bigfrom;
    double pred_l, sim_l, adjrate;
    TransMat diag, Spopmat, Rpopmat;

    if(by_population != 0)
    {
	if(by_population == 1) 
	{
	    diag.SetSize(num_stages_);
	    Spopmat.SetSize(num_stages_);
	    Rpopmat.SetSize(num_stages_);

	    for(i = 0; i < num_habitats_; i++)
	    {
		for(j = 0; j < num_stages_; j++)
		{
		    for(k = 0; k < num_stages_; k++)
		    {
			bigto = (i * num_stages_) + k;
			bigfrom = (i * num_stages_) + j;

			Spopmat.SetElement(k, j, host_S_[num_classes_ * bigto + bigfrom]);
			Rpopmat.SetElement(k, j, host_R_[num_classes_ * bigto +  bigfrom]);
		    }
		}

		pred_l = (Spopmat + Rpopmat).Lambda();
		sim_l = (Spopmat * (Rpopmat + diag)).Lambda();
		adjrate = pred_l / sim_l;

		for(l = (i * num_stages_); l < ((i * num_stages_) + num_stages_); l++)
		{
		    CarryState(int(round(double(CountIndividualsInClass(l)) * adjrate)), l);
		}
	    }
	}
	else
	{
	    assert(1 == 0);
/*
	    diag.SetSize(num_classes_);
	    pred_l = (S[e]+R[e]).Lambda();
	    sim_l = (S[e]*(R[e]+diag)).Lambda();
	    adjrate = pred_l/sim_l;
	    
	    for (i=0;i<(s*nhab);i++)
	    {
		CarryState(int(round(double(I[i].size())*adjrate)),i);
	    }
*/
	} // if(by_population == 1)
    } // if(by_population != 0)

    timer_.EndMethod();
}

void Landscape_gpu::HabCarry(const int k)
{
    std::vector<double> prop(num_habitats_);

    for(int h = 0; h < num_habitats_; h++)
    {
	if (k < 0)
	{
	    prop[h] = double(habitat_carrying_capacities_[h]) / double(CountIndividualsInHabitat(h));
	}
	else
	{
	    prop[h] = double(k) / double(CountIndividualsInHabitat(h));
	}

	if (prop[h] > 1) 
	{
	    prop[h] = 1.0;
	}
    }

    for(int j = 0; j < num_classes_; j++)
    {
	int habitat = j / num_stages_;
	CarryState(prop[habitat] * CountIndividualsInClass(j), j);
    }
}

void Landscape_gpu::LandCarry()
{
    const int pop_size = CountIndividuals();
    const double pr = static_cast<double>(carrying_capacity_) / pop_size;

    for (int j = 0; j < num_classes_; j++)
    {
	CarryState(pr * CountIndividualsInClass(j), j);
    }
}

void Landscape_gpu::CarryState(const int max_size, const int state)
{
    timer_.BeginMethod("CarryState");

    const int total_pop_size = CountIndividuals();
    const int state_count = CountIndividualsInClass(state);

    if (max_size < state_count)
    {
	//XXX:remove this line, just for testing
//	CompactIndividuals();

	int num_to_delete = state_count - max_size;

	thrust::fill(shuffling_dartboard_.begin(), shuffling_dartboard_.end(),
	    -1);

	kernel_params_.state = state;
	kernel_params_.num_individuals = total_pop_size;

	cuda::FillShuffleVector(kernel_params_);

	const int state_size_check = thrust::count_if(
	    shuffling_dartboard_.begin(), shuffling_dartboard_.end(), 
	    functors::non_negative_int());

	assert(state_size_check == state_count);

	thrust::copy_if(shuffling_dartboard_.begin(),
	    shuffling_dartboard_.end(), shuffling_dartboard_.begin(), 
	    functors::non_negative());

	thrust::fill(shuffling_dartboard_.begin() + num_to_delete,
	    shuffling_dartboard_.end(), -1);

	thrust::fill(
	    thrust::make_permutation_iterator(
		classes_.begin(), 
		shuffling_dartboard_.begin()),
	    thrust::make_permutation_iterator(
		classes_.begin() + num_to_delete, 
		shuffling_dartboard_.begin() + num_to_delete), -1);

	//const int pop_size_after = CountIndividualsInClass(state);
	//assert(pop_size_after == max_size);

	CompactIndividuals();
    }

    timer_.EndMethod();
}

void Landscape_gpu::Advance()
{
    current_generation_++;
}

// Initialization

void Landscape_gpu::CalculateCapacities()
{
    carrying_capacity_ = 0;

    for(int i = 0; i < num_habitats_; i++)
    {
	carrying_capacity_ += habitat_carrying_capacities_[i];
    }

    max_size_ = carrying_capacity_ * LANDSCAPE_VECTOR_SIZE_MULTIPLIER;
}

void Landscape_gpu::InitializeKernelParameters()
{
    // Landscape parameters
    kernel_params_.max_individuals = max_size_;
    kernel_params_.num_classes = num_classes_;
    kernel_params_.num_habitats = num_habitats_;
    kernel_params_.num_stages = num_stages_;
    kernel_params_.selfing_rate = selfing_rate_;
    kernel_params_.multiple_paternity = multiple_paternity_;

    // Individual vectors
    kernel_params_.generations_born = thrust::raw_pointer_cast(
	generations_born_.data());
    kernel_params_.generations_last_changed = thrust::raw_pointer_cast(
	generations_last_changed_.data());
    kernel_params_.generations_last_reproduced = thrust::raw_pointer_cast(
	generations_last_reproduced_.data());
    kernel_params_.classes = thrust::raw_pointer_cast(classes_.data());
    kernel_params_.ids = thrust::raw_pointer_cast(ids_.data());
    kernel_params_.maternal_ids = thrust::raw_pointer_cast(
	maternal_ids_.data());
    kernel_params_.paternal_ids = thrust::raw_pointer_cast(
	paternal_ids_.data());
    kernel_params_.num_offspring = thrust::raw_pointer_cast(
	num_offspring_.data());
    kernel_params_.genotypes = thrust::raw_pointer_cast(
	genotypes_.data());

    // Shuffling vectors
    kernel_params_.shuffling_dartboard = thrust::raw_pointer_cast(
	shuffling_dartboard_.data());

    // Offspring vectors
    kernel_params_.offspring_maternal_indices = thrust::raw_pointer_cast(
	offspring_maternal_indices_.data());
    kernel_params_.offspring_paternal_indices = thrust::raw_pointer_cast(
	offspring_paternal_indices_.data());
    kernel_params_.offspring_classes = thrust::raw_pointer_cast(
	offspring_classes_.data());
    kernel_params_.offspring_genotypes = thrust::raw_pointer_cast(
	offspring_genotypes_.data());

    // Transition matrices
    kernel_params_.S = thrust::raw_pointer_cast(S_.data());
    kernel_params_.R = thrust::raw_pointer_cast(R_.data());
    kernel_params_.M = thrust::raw_pointer_cast(M_.data());

    // RNG state vector
    kernel_params_.rand_states = thrust::raw_pointer_cast(
	rand_states_.data());
}

void Landscape_gpu::InitializeRandStates()
{
    PRNG rng(rand_seed_);
    thrust::host_vector<unsigned long> h_per_thread_seeds(max_size_);
	
    for(int i = 0; i < max_size_; i++)
    {
	h_per_thread_seeds[i] = rng.RandomUInt();
    }

    thrust::device_vector<unsigned long> d_per_thread_seed = h_per_thread_seeds;
    cuda::InitializeRandStates(thrust::raw_pointer_cast(d_per_thread_seed.data()), kernel_params_);
}

void Landscape_gpu::ResizeVectors()
{   
    generations_born_.resize(max_size_);
    generations_last_changed_.resize(max_size_);
    generations_last_reproduced_.resize(max_size_);
    classes_.resize(max_size_);
    ids_.resize(max_size_);
    maternal_ids_.resize(max_size_);
    paternal_ids_.resize(max_size_);
    num_offspring_.resize(max_size_);
    genotypes_.resize(max_size_);

    shuffling_dartboard_.resize(max_size_);

    offspring_maternal_indices_.resize(OFFSPRING_ALLOC_CHUNK_SIZE);
    offspring_paternal_indices_.resize(OFFSPRING_ALLOC_CHUNK_SIZE);
    offspring_classes_.resize(OFFSPRING_ALLOC_CHUNK_SIZE);
    offspring_genotypes_.resize(OFFSPRING_ALLOC_CHUNK_SIZE);

    rand_states_.resize(max_size_);
}

// Landscape_gpu<=>Landscape Individuals Transfer

void Landscape_gpu::ExtractIndividualsFromLandscape(Landscape_statistics &landscape, std::vector<CudaIndividual> &individuals)
{
    for(int demo_class = 0; demo_class < num_classes_; demo_class++)
    {
	landscape.resetStage(demo_class);
	PackedIndividual ind = landscape.getNextInd(demo_class);

	while(ind.cl != -1)
	{
	    landscape.advanceStagePtr(demo_class);

	    CudaIndividual cuda_ind = { ind.id, ind.mid, ind.pid, ind.cl, ind.gen, ind.changed, ind.lastrep };
	    individuals.push_back(cuda_ind);
	    
	    ind = landscape.getNextInd(demo_class);
	}
    }
}

void Landscape_gpu::InjectIndividualsIntoLandscape(const std::list<CudaIndividual> &individuals, Landscape_statistics &landscape)
{
    std::vector<int> pop_sizes(num_classes_);

    for(int i = 0; i < num_classes_; i++)
    {
	pop_sizes[i] = CountIndividualsInClass(i);
    }

    landscape.popsizeset(pop_sizes);

    landscape.GCAlleles();

    StepAlleleTbl *s = new StepAlleleTbl;
    landscape.Atbl_push_back(s);

    std::list<CudaIndividual>::const_iterator iter = individuals.begin();

    while(iter != individuals.end())
    {
	const CudaIndividual &cuda_ind = *iter;

	assert(cuda_ind.cl >= 0);
	assert(cuda_ind.cl < num_classes_);

	PackedIndividual ind;

	ind.id = cuda_ind.id;
	ind.mid = cuda_ind.mid;
	ind.pid = cuda_ind.pid;
	ind.cl = cuda_ind.cl;
	ind.gen = cuda_ind.gen;
	ind.changed = cuda_ind.changed;
	ind.lastrep = cuda_ind.lastrep;

	std::cout << ind << std::endl;

	landscape.addIndividual(ind, -1);

	iter++;
    }
}

// Getters/Setters

void Landscape_gpu::SetIndividuals(const std::vector<CudaIndividual> &h_individuals)
{
//    printf("SetIndividuals():\n");

    thrust::fill(classes_.begin(), classes_.end(), -1);
    int num_to_transfer = h_individuals.size();
    int num_blocks = 0;

//    printf("Transferring %d individuals to GPU\n", num_to_transfer);

    if(num_to_transfer > MAX_INDIVIDUAL_TRANSFER_SIZE)
    {
//	printf("Large number of individuals to transfer.\n");
//	printf("Partitioning individuals into smaller blocks of size %d.\n",
//	    MAX_INDIVIDUAL_TRANSFER_SIZE);

	thrust::device_vector<CudaIndividual> d_individuals(MAX_INDIVIDUAL_TRANSFER_SIZE);
	int num_transferred = 0;

	while(num_transferred < num_to_transfer)
	{
	    num_blocks++;

	    int transfer_size = num_to_transfer > MAX_INDIVIDUAL_TRANSFER_SIZE ? MAX_INDIVIDUAL_TRANSFER_SIZE : num_to_transfer;
	    d_individuals.assign(h_individuals.begin() + num_transferred, h_individuals.begin() + num_transferred + transfer_size);

//	    printf("Block #%d: size=%d\n", num_blocks, transfer_size);
	    
	    thrust::transform(
		d_individuals.begin(), 
		d_individuals.begin() + transfer_size,
		thrust::make_zip_iterator(
		    thrust::make_tuple(
			ids_.begin() + num_transferred,
			maternal_ids_.begin() + num_transferred,
			paternal_ids_.begin() + num_transferred,
			classes_.begin() + num_transferred,
			generations_born_.begin() + num_transferred,
			generations_last_changed_.begin() + num_transferred,
			generations_last_reproduced_.begin() + num_transferred)),
		functors::individidual_to_tuple());
	}
    }
    else
    {
	thrust::device_vector<CudaIndividual> d_individuals(h_individuals.begin(), h_individuals.begin() + num_to_transfer);
	
	thrust::transform(
	    d_individuals.begin(), 
	    d_individuals.end(), 
	    thrust::make_zip_iterator(
		thrust::make_tuple(
		    ids_.begin(),
		    maternal_ids_.begin(),
		    paternal_ids_.begin(),
		    classes_.begin(),
		    generations_born_.begin(),
		    generations_last_changed_.begin(),
		    generations_last_reproduced_.begin())),
	    functors::individidual_to_tuple());
    }

    CompactIndividuals();
}

void Landscape_gpu::GetIndividuals(std::list<CudaIndividual> &h_individuals)
{
    const int pop_size = CountIndividuals();

    thrust::device_vector<CudaIndividual> d_individuals(pop_size);

    thrust::transform(
	thrust::make_zip_iterator(
	    thrust::make_tuple(
		ids_.begin(),
		maternal_ids_.begin(),
		paternal_ids_.begin(),
		classes_.begin(),
		generations_born_.begin(),
		generations_last_changed_.begin(),
		generations_last_reproduced_.begin())),
	thrust::make_zip_iterator(
	    thrust::make_tuple(
		ids_.begin() + pop_size,
		maternal_ids_.begin() + pop_size,
		paternal_ids_.begin() + pop_size,
		classes_.begin() + pop_size,
		generations_born_.begin() + pop_size,
		generations_last_changed_.begin() + pop_size,
		generations_last_reproduced_.begin() + pop_size)),
	d_individuals.begin(),
	functors::tuple_to_individual());

    thrust::host_vector<CudaIndividual> th_individuals(pop_size);
    thrust::copy(d_individuals.begin(), d_individuals.end(), th_individuals.begin());

    h_individuals.assign(th_individuals.begin(), th_individuals.end());
}

// Private Simulation Methods

int Landscape_gpu::CountIndividuals()
{
    return thrust::count_if(classes_.begin(), classes_.end(), functors::non_negative());
}

int Landscape_gpu::CountIndividualsInHabitat(int habitat)
{
    assert(habitat >= 0);
    assert(habitat < num_habitats_);

    int min_class = habitat * num_stages_;
    int max_class = min_class + num_stages_;

    return thrust::count_if(classes_.begin(), classes_.end(), functors::in_range(min_class, max_class));
}

int Landscape_gpu::CountIndividualsInClass(int cl)
{
    assert(cl >= 0);
    assert(cl < num_classes_);

    return thrust::count(classes_.begin(), classes_.end(), cl);
}

void Landscape_gpu::CalculateMaleGameteMatrix()
{
    std::vector<double> n(num_classes_ * num_classes_);
    std::vector<int> class_sizes(num_classes_);

    thrust::host_vector<double> M(num_classes_ * num_classes_);

    for(int cl = 0; cl < num_classes_; cl++)
    {
	class_sizes[cl] = CountIndividualsInClass(cl);
    }

    for(int to = 0; to < num_classes_; to++)
    {
	int column_start_index = to * num_classes_;
	double sum_weighted_class_size = 0.;

	for(int from = 0; from < num_classes_; from++)
	{
	    int linear_index = column_start_index + from;

	    int m_value = host_M_[linear_index];
	    int class_size = double(class_sizes[from]);

	    n[linear_index] = (m_value * class_size);
	    sum_weighted_class_size += n[linear_index];
	}

	if(sum_weighted_class_size > 0)
	{
	    double sum_probability = 0.;

	    for(int from = 0; from < num_classes_; from++)
	    {
		int linear_index = column_start_index + from;

		double probability = n[linear_index] / sum_weighted_class_size;

		M[linear_index] = probability;
		sum_probability =+ probability;
	    }

	    if (sum_probability > 1.)
	    {
		if (sum_probability > 1.1) //something is very wacky and the program should terminate
		{
		    cerr << "The probabilities of choosing a male gamete class total to more than 1: total = "<< sum_probability << endl;
		    assert(sum_probability <= 1);
		}
		else
		{
		    for(int from = 0; from < num_classes_; from++)
		    {
			int linear_index = column_start_index + from;
			M[linear_index] = M[linear_index] / sum_probability;
		    }
		}
	    }
	}
    }

    M_ = M;
}

void Landscape_gpu::SetTransitionMatrices()
{
    S_ = host_S_;
    R_ = host_R_;
}

void Landscape_gpu::CompactIndividuals()
{
    const int pop_size_before = CountIndividuals();

    thrust::copy_if(
	thrust::make_zip_iterator(
	    thrust::make_tuple(
		ids_.begin(),
		maternal_ids_.begin(),
		paternal_ids_.begin(),
		classes_.begin(),
		generations_born_.begin(),
		generations_last_changed_.begin(),
		generations_last_reproduced_.begin(),
		genotypes_.begin())), 
	thrust::make_zip_iterator(
	    thrust::make_tuple(
		ids_.end(),
		maternal_ids_.end(),
		paternal_ids_.end(),
		classes_.end(),
		generations_born_.end(),
		generations_last_changed_.end(),
		generations_last_reproduced_.end(),
		genotypes_.end())),
	classes_.begin(),
	thrust::make_zip_iterator(
	    thrust::make_tuple(
		ids_.begin(),
		maternal_ids_.begin(),
		paternal_ids_.begin(),
		classes_.begin(),
		generations_born_.begin(),
		generations_last_changed_.begin(),
		generations_last_reproduced_.begin(),
		genotypes_.begin())),
	functors::non_negative());

    thrust::fill(classes_.begin() + pop_size_before, classes_.end(), -1);

    assert(pop_size_before == CountIndividuals());
}

// Debug Methods

void Landscape_gpu::PrintSimulationStatistics()
{
    timer_.PrintStatistics();
}

void Landscape_gpu::PrintIndividuals(const int num_to_print)
{
    kernel_params_.num_individuals = (num_to_print < max_size_) ? num_to_print : max_size_;
    cuda::PrintIndividuals(kernel_params_);
}

void Landscape_gpu::PrintLandscape()
{
    PrintIndividuals(max_size_);
}

} /* namespace rmetasim_gpu */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

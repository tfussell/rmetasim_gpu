
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
#include <time.h>
#include <stdint.h>
#ifdef __Linux__
#include <sys/time.h>
#else
#include <Windows.h>
#endif
#include <cmath>
#include <stdio.h>

/* Local Includes */
#include "gpu/SimulationTimer.h"

namespace rmetasim_gpu {

double Mean(const std::vector<double> &vec)
{
    double sum = 0;

    for(size_t i = 0; i < vec.size(); i++)
    {
	sum += vec[i];
    }

    double mean = sum / vec.size();

    return mean;
}

double StandardDeviation(const std::vector<double> &vec)
{
    double mean = Mean(vec);
    double sum_squared_differences = 0;

    for(size_t i = 0; i < vec.size(); i++)
    {
	double diff = mean - vec[i];
	sum_squared_differences += diff * diff;
    }

    double variance = sum_squared_differences / vec.size();

    return std::sqrt(variance);
}

double Min(const std::vector<double> &vec)
{
    double min = vec[0];

    for(size_t i = 0; i < vec.size(); i++)
    {
	if(vec[i] < min)
	{
	    min = vec[i];
	}
    }

    return min;
}

double Max(const std::vector<double> &vec)
{
    double max = vec[0];

    for(size_t i = 0; i < vec.size(); i++)
    {
	if(vec[i] > max)
	{
	    max = vec[i];
	}
    }

    return max;
}

SimulationTimer::SimulationTimer() : method_name_(""), method_start_time_(-1)
{
#ifdef __Linux__
    zeroClock = clock();
    gettimeofday(&start, NULL);
#else
    LARGE_INTEGER ticksPerSecond;
    QueryPerformanceFrequency(&ticksPerSecond);
    updateFrequency = double(ticksPerSecond.QuadPart) / 1000000.0;
#endif
}

SimulationTimer::~SimulationTimer()
{

}

void SimulationTimer::PrintStatistics()
{
    TimingMap::iterator iter;

    for(iter = method_timings_.begin(); iter != method_timings_.end(); iter++)
    {
	const std::string method = iter->first;
	const std::vector<double> &timings = iter->second;

	if(timings.size() > 0)
	{
	    if(timings.size() == 1)
	    {
		printf("%s: %.3fus\n", method.c_str(), timings[0]);
	    }
	    else
	    {
#if 0
		printf("%s:\n", method.c_str());
		printf("\tMin: %.3fus\n", Min(timings));
		printf("\tMax: %.3fus\n", Max(timings));
		printf("\tAverage: %.3fus\n", Mean(timings));
		printf("\tStandard Deviation: %.3fus\n", StandardDeviation(timings));
#else
		printf("%s: %.3fus\n", method.c_str(), Mean(timings));
#endif
	    }
	}
    }
}

void SimulationTimer::BeginMethod(const std::string &method_name)
{
    assert(method_name != "");
    assert(method_start_time_ == -1);
    assert(method_name_ == "");

    TimingMap::iterator iter = method_timings_.find(method_name);

    if(iter == method_timings_.end())
    {
	method_timings_[method_name] = std::vector<double>();
    }

    method_name_ = method_name;
    method_start_time_ = GetTime();
}

void SimulationTimer::EndMethod()
{
    assert(method_start_time_ >= 0);
    assert(method_name_ != "");

    double end_time = GetTime();
    double elapsed_time = end_time - method_start_time_;

    method_timings_[method_name_].push_back(elapsed_time);

    method_start_time_ = -1;
    method_name_ = "";
}

double SimulationTimer::GetTime()
{
#ifdef __Linux__
    struct timeval now;
    gettimeofday(&now, NULL);
    return (now.tv_sec - start.tv_sec) * 1000000 + (now.tv_usec - start.tv_usec);
#else
    LARGE_INTEGER tick;
    QueryPerformanceCounter(&tick);
    return static_cast<unsigned long>((tick.QuadPart - start) / updateFrequency);
#endif
}

} /* namespace rmetasim_gpu */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

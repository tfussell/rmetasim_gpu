
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

#ifndef __SimulationTimer_H__
#define __SimulationTimer_H__

/* External Includes */
#include <map>
#include <vector>
#include <string>
#include <stdint.h>

namespace rmetasim_gpu {

/* Typedefs */
typedef std::map<std::string, std::vector<double> > TimingMap;

/* Object Declarations */
class SimulationTimer
{
public:
    SimulationTimer();
    ~SimulationTimer();

    void PrintStatistics();

    void BeginMethod(const std::string &method);
    void EndMethod();

private:
    double GetTime();

    std::string method_name_;
    double method_start_time_;

    TimingMap method_timings_;

#ifdef __Linux__
    struct timeval start;
    clock_t zeroClock;
#else /* windows */
    int64_t start;
    double updateFrequency;
#endif
}; /* class SimulationTimer */

} /* namespace rmetasim_gpu */

#endif /* #ifndef __SimulationTimer_H__ */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

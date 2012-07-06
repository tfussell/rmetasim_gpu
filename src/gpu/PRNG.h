
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

#ifndef __PRNG_H__
#define __PRNG_H__

/* External Includes */
#include <iostream>

namespace rmetasim_gpu {

/* Object Declarations */
class PRNG
{
public:
    PRNG(const int seed);
    PRNG();
    ~PRNG();

    void SeedRNG(const int seed);

    int RandomInt(int min, int max);
    int RandomInt(int max);

    double RandomDouble();
    double RandomDouble(double min);
    double RandomDouble(double min, double max);

    unsigned int RandomUInt();

    friend std::ostream &operator<<(std::ostream &out, const PRNG &rng);
    friend std::istream &operator>>(std::istream &out, PRNG &rng);

private:
    unsigned int x[5];
}; /* class PRNG */

} /* namespace rmetasim_gpu */

#endif /* ifndef __PRNG_H__ */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/


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
#include "gpu/PRNG.h"

namespace rmetasim_gpu {

PRNG::PRNG()
{
    SeedRNG(0);
}

PRNG::PRNG(const int seed)
{
    SeedRNG(seed);
}

PRNG::~PRNG()
{

}

void PRNG::SeedRNG(const int seed)
{
    int i;
    unsigned int s = seed;

    for(i = 0; i < 5; i++)
    {
	s = s * 29943829 - 1;
	x[i] = s;
    }

    for(i = 0; i < 19; i++)
    {
	RandomUInt();
    }
}

unsigned int PRNG::RandomUInt()
{
    unsigned long long int sum = 0L;

    sum = (unsigned long long int)2111111111UL * (unsigned long long int)(x[3]) +
          (unsigned long long int)1492 * (unsigned long long int)(x[2]) +
  	  (unsigned long long int)1776 * (unsigned long long int)(x[1]) +
	  (unsigned long long int)5115 * (unsigned long long int)(x[0]) +
	  (unsigned long long int)(x[4]);

    x[3] = x[2];
    x[2] = x[1];
    x[1] = x[0];
    x[4] = (unsigned int)(sum >> 32);
    x[0] = (unsigned int)sum;

    return x[0];
}

double PRNG::RandomDouble()
{
    return (double) RandomUInt() * (1. / (65536. * 65536.));
}

double PRNG::RandomDouble(double max)
{
    return RandomDouble() * max;
}

double PRNG::RandomDouble(double min, double max)
{
    return min + RandomDouble() * (max - min);
}

int PRNG::RandomInt(int max)
{
    return RandomInt(0, max);
}

int PRNG::RandomInt(int min, int max)
{
    if(max <= min)
    {
	if(max == min)
        {
	    return min;
        }
	else
        {
	    return 0x80000000;
        }
    }

    unsigned int interval = (unsigned int)(max - min);
    return ((int)(interval * RandomDouble())) + min;
}

std::istream &operator>>(std::istream &in, PRNG &rng)
{
    for(int i = 0; i < 5; i++)
    {
	in >> rng.x[i];
    }

    return in;
}

std::ostream &operator<<(std::ostream &out, const PRNG &rng)
{
    out << rng.x[0] << " " << rng.x[1] << " " << rng.x[2] << " " << rng.x[3] << " " << rng.x[4];

    return out;
}

} /* namespace rmetasim_gpu */

/*
  ;;; Local Variables:        ***
  ;;; mode: C++               ***
  ;;; minor-mode:  font-lock  ***
  ;;; End:                    ***
*/

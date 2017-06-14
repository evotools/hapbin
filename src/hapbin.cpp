/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014-2017 Colin MacLean <cmaclean@illinois.edu>
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
 *
 */

#include "hapbin.hpp"
#include <vector>
#include <sstream>
#include <unistd.h>

double binom_2(double n)
{
    return n*(n-1.0)*0.5;
}

double nearest(double target, double number)
{
    if (target == 0)
        return number;
    return std::round(number/target)*target;
}

Stats stats(const std::vector<double>& list)
{
    Stats s;
    if (list.size() == 0)
        return s;

    double total = 0.0;
    for (double v : list)
    {
        total += v;
    }

    s.mean = total/(double)list.size();

    double sqtotal = 0.0;
    for (double v : list)
    {
        sqtotal += std::pow(v - s.mean, 2);
    }

    s.stddev = std::sqrt(sqtotal/(double)list.size());

    return s;
}

std::vector<std::string> splitString(const std::string input, char delim)
{
    std::stringstream ss(input);
    std::vector<std::string> split;
    std::string part;
    while (std::getline(ss,part,delim))
        split.push_back(part);
    return split;
}

void filter(unsigned long long* in, unsigned long long* out, std::size_t inlen, std::size_t outlen, const PopKey& pk)
{
    assert(outlen <= inlen);
    assert(out[0] == 0ULL);
    std::size_t outpos = 0ULL;
    constexpr std::size_t bits = sizeof(unsigned long long)*8;
    for(size_t i = 0; i < inlen; ++i)
    {
        bool state = (in[i/bits] & (1ULL << (i % bits))) > 0;
        if (pk[i])
        {
            if (state)
                out[outpos/bits] |= (1ULL << (outpos % bits));
            ++outpos;
        }
    }
}


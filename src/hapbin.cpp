/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014  Colin MacLean <s0838159@sms.ed.ac.uk>
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

#if defined(__MINGW32__) && !defined(_ISOC11_SOURCE)
void* aligned_alloc(size_t alignment, size_t size)
{
    return __mingw_aligned_malloc(size, alignment);
}
void aligned_free(void* ptr)
{
    __mingw_aligned_free(ptr);
}
#elif (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !defined(_ISOC11_SOURCE)
void* aligned_alloc(size_t alignment, size_t size)
{
    void* ret;
    if (posix_memalign(&ret, alignment, size) != 0)
        throw std::bad_alloc();
    return ret;
}
#endif

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

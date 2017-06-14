/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2015-2017 Colin MacLean <cmaclean@illinois.edu>
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

#ifndef IHSPAIRFINDER_H
#define IHSPAIRFINDER_H

#include "ehhpairfinder.hpp"
#include <vector>
#include <mutex>

class IhsPairJob;
class IhsPairFinder
{
public:
    IhsPairFinder(const HapMap* hm, double cutoff, double minMaf, double scale, unsigned long long window, std::size_t maxBreadth);
    IhsPairJob* calcRange(std::size_t start, std::size_t end);
    void add(const std::vector<EHHPair>& v);
    std::size_t numPairs();
    std::vector<EHHPair> results() const { return m_results; }
protected:
    const HapMap* m_hm;
    const double m_cutoff;
    const double m_minMaf;
    const double m_scale;
    const unsigned long long m_window;
    std::size_t  m_maxBreadth;
    std::vector<EHHPair> m_results;
    std::mutex m_mutex;
};

#endif // IHSPAIRFINDER_H

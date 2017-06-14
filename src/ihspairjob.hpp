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

#ifndef IHSPAIRJOB_H
#define IHSPAIRJOB_H

#include "ehh.hpp"
#include "hapmap.hpp"

class IhsPairJob
{
public:
    IhsPairJob();
    IhsPairJob(std::size_t start, std::size_t end);
    void add(const EHHPair &result);
    int load(const std::string& file);
    void save(const std::string& file);
    std::vector<EHHPair> results() const;
    void saveAscii(const std::string& out, const HapMap& hm);
    std::size_t size() const { return m_results.size(); }
    std::size_t start() const { return m_range[0]; }
    std::size_t end() const { return m_range[1]; }
protected:
    std::size_t m_range[2];
    double m_sumIhs[4];
    std::vector<EHHPair> m_results;
    static const uint64_t magicNumber;
};

#endif // IHSPAIRJOB_H

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

#include "ihspairfinder.hpp"
#include "ehh.hpp"
#include "ihspairjob.hpp"

IhsPairFinder::IhsPairFinder(const HapMap* hm, double cutoff, double minMaf, double scale, unsigned long long window, std::size_t maxBreadth)
    : m_hm(hm)
    , m_cutoff(cutoff)
    , m_minMaf(minMaf)
    , m_scale(scale)
    , m_window(window)
    , m_maxBreadth(maxBreadth)
{}

std::size_t IhsPairFinder::numPairs()
{
    std::size_t ret = 0;
    EhhPairFinder pf(m_hm, m_cutoff, m_minMaf, m_scale, m_maxBreadth);
    for(std::size_t focus1 = 0; focus1 < m_hm->numSnps(); ++focus1)
    {
        //if (focus1 % 100 == 0)
        std::cout << focus1 << " " << ret << std::endl;
        unsigned long long focus1pos = m_hm->physicalPosition(focus1);
        std::size_t focus2 = focus1+1;
        if (focus2 >= m_hm->numSnps())
            continue;
        unsigned long long focus2pos = m_hm->physicalPosition(focus2);
        while (focus2pos-focus1pos < m_window)
        {
            if (pf.passesAf(focus1, focus2))
                ++ret;
            ++focus2;
            if (focus2 >= m_hm->numSnps())
                break;
            focus2pos = m_hm->physicalPosition(focus2);
        }
    }
    return ret;
}

IhsPairJob* IhsPairFinder::calcRange(std::size_t start, std::size_t end)
{
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    IhsPairJob *job = new IhsPairJob(start, end);
    static unsigned long long total_pairs = 0.0;
    #pragma omp parallel shared(start, end, pairs)
    {
        EhhPairFinder pf(m_hm, m_cutoff, m_minMaf, m_scale, m_maxBreadth);
        #pragma omp single
        {
            for(std::size_t focus1 = start; focus1 < end; ++focus1)
            {
                unsigned long long focus1pos = m_hm->physicalPosition(focus1);
                std::size_t focus2 = focus1+1;
                if (focus2 >= m_hm->numSnps())
                    continue;
                unsigned long long focus2pos = m_hm->physicalPosition(focus2);
                while (focus2pos-focus1pos < m_window)
                {
                    if (pf.passesAf(focus1, focus2))
                        pairs.push_back({focus1,focus2});
                    ++focus2;
                    if (focus2 >= m_hm->numSnps())
                        break;
                    focus2pos = m_hm->physicalPosition(focus2);
                }
            }
            std::cout << "Calculating " << pairs.size() << " pairs for range " << start << "-" << end << std::endl;
        } // implicit barrier
        #pragma omp for schedule(dynamic,1)
        for(auto i = pairs.begin(); i < pairs.end(); ++i)
        {
            auto res = pf.calcEhhPair(i->first, i->second, false);
            if (res.focus[0] != 0UL && res.focus[1] != 0UL)
            {
                m_mutex.lock();
                //m_results.push_back(res);
                job->add(res);
                m_mutex.unlock();
            }
        }
        #pragma omp single
        {
            total_pairs += job->size();
            std::cout << "Calculated: " << job->size() << ". Total pairs: " << total_pairs << std::endl;
        }
    }
    return job;
}

void IhsPairFinder::add(const std::vector<EHHPair>& v)
{
    for (const auto& i : v)
        m_results.push_back(i);
}

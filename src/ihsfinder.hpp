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

#ifndef IHSFINDER_H
#define IHSFINDER_H
#include "ehhfinder.hpp"
#include <map>
#include <mutex>
#ifdef __MINGW32__
//from https://github.com/meganz/mingw-std-threads
#include <windows.h>
#include "mingw.mutex.h"
#endif
#include <functional>
#include <thread>
#include <atomic>

class IHSFinder
{
public:
    using IndexMap = std::map<std::size_t, double>;
    using IhsInfoMap = std::map<std::size_t, IhsScore>;
    using XIhsInfoMap = std::map<std::size_t, XPEHH>;
    using FreqVecMap = std::map<double, std::vector<double>>;
    using StatsMap = std::map<double, Stats>;
    IHSFinder(std::size_t snpLength, double cutoff, double minMAF, double scale, unsigned long long maxExtend, int bins);

    FreqVecMap unStdIHSByFreq() const { return m_unStdIHSByFreq; }
    FreqVecMap unStdNSLByFreq() const { return m_unStdNSLByFreq; }
    IhsInfoMap unStdByIndex() const { return m_unStdByIndex; }
    XIhsInfoMap unStdXIHSByIndex() const { return m_unStdXIHSByIndex; }
    IndexMap    freqsByIndex() const    { return m_freqsByIndex; }
    IhsInfoMap    std() const { return m_std; }

    unsigned long long numCompleted() const { return m_counter; }
    unsigned long long numReachedEnd() const { return m_reachedEnd; }
    unsigned long long numOutsideMaf() const { return m_outsideMaf; }
    unsigned long long numNanResults() const { return m_nanResults; }

    template <bool Binom>
    void run(const HapMap* map, std::size_t start, std::size_t end);
    template <bool Binom>
    void runXpehh(const HapMap* mA, const HapMap* mB, std::size_t start, std::size_t end);
    void addData(const IHSFinder::IndexMap& freqsByIndex,
                 const IHSFinder::IhsInfoMap& unStandIHSByIndex,
                 const IHSFinder::FreqVecMap& unStdIHSByFreq,
                 const IHSFinder::FreqVecMap& unStdNSLByFreq,
                 long long unsigned int reachedEnd,
                 long long unsigned int outsideMaf,
                 long long unsigned int nanResults);
    void addXData(const IndexMap& freqsBySite,
                  const XIhsInfoMap& unStandXIHSByIndex,
                  const FreqVecMap& unStandIHSByFreq,
                  unsigned long long reachedEnd,
                  unsigned long long outsideMaf,
                  unsigned long long nanResults);
    void normalize();

protected:
    void processEHH(const EHH& ehh, std::size_t index);
    void processXPEHH(XPEHH& e, size_t index);
    std::size_t m_snpLength;
    double m_cutoff;
    double m_minMAF;
    double m_scale;
    unsigned long long m_maxExtend;
    int m_bins;
    std::mutex m_mutex;

    IndexMap    m_freqsByIndex;
    IhsInfoMap  m_unStdByIndex;
    XIhsInfoMap m_unStdXIHSByIndex;
    FreqVecMap  m_unStdIHSByFreq;
    FreqVecMap  m_unStdNSLByFreq;
    IhsInfoMap  m_std;
    std::atomic<unsigned long long> m_counter;
    std::atomic<unsigned long long> m_reachedEnd;
    std::atomic<unsigned long long> m_outsideMaf;
    std::atomic<unsigned long long> m_nanResults;
};

#include "ihsfinder-impl.hpp"

#endif // IHSFINDER_H

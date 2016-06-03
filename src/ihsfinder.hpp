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
    using LineMap = std::map<std::size_t, double>;
    using IhsInfoMap = std::map<std::size_t, IhsScore>;
    using XIhsInfoMap = std::map<std::size_t, XPEHH>;
    using FreqVecMap = std::map<double, std::vector<double>>;
    using StatsMap = std::map<double, Stats>;
    
    IHSFinder(std::size_t snpLength, double cutoff, double minMAF, double scale, unsigned long long maxExtend, int bins);
    FreqVecMap unStdIHSByFreq() const { return m_unStandIHSByFreq; }
    IhsInfoMap unStdIHSByLine() const { return m_unStandIHSByLine; }
    XIhsInfoMap unStdXIHSByLine() const { return m_unStandXIHSByLine; }
    LineMap    freqsByLine() const    { return m_freqsByLine; }
    unsigned long long numCompleted() const { return m_counter; }
    unsigned long long numReachedEnd() const { return m_reachedEnd; }
    unsigned long long numOutsideMaf() const { return m_outsideMaf; }
    unsigned long long numNanResults() const { return m_nanResults; }
    
    void run(HapMap* map, std::size_t start, std::size_t end);
    void runXpehh(HapMap* mA, HapMap* mB, std::size_t start, std::size_t end);
    LineMap normalize();
    
    void addData(const LineMap& freqsBySite, const IhsInfoMap& unStandIHSByLine, const FreqVecMap& unStandIHSByFreq, unsigned long long reachedEnd, unsigned long long outsideMaf, unsigned long long nanResults);
    void addXData(const LineMap& freqsBySite, const XIhsInfoMap& unStandXIHSByLine, const FreqVecMap& unStandIHSByFreq, unsigned long long reachedEnd, unsigned long long outsideMaf, unsigned long long nanResults);
    
protected:
    void processEHH(const EHH& ehh, std::size_t line);
    void processXPEHH(XPEHH& e, size_t line);
    
    std::size_t m_snpLength;
    double m_cutoff;
    double m_minMAF;
    double m_scale;
    unsigned long long m_maxExtend;
    int m_bins;
    
    std::mutex m_mutex;
    std::mutex m_freqmutex;
    LineMap    m_freqsByLine;
    IhsInfoMap m_unStandIHSByLine;
    XIhsInfoMap m_unStandXIHSByLine;
    FreqVecMap m_unStandIHSByFreq;
    LineMap    m_standIHSSingle;
    
    std::atomic<unsigned long long> m_counter;
    std::atomic<unsigned long long> m_reachedEnd;
    std::atomic<unsigned long long> m_outsideMaf;
    std::atomic<unsigned long long> m_nanResults;
};


#endif // IHSFINDER_H

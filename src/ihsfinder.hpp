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
    using FreqVecMap = std::map<double, std::vector<double>>;
    using StatsMap = std::map<double, Stats>;
    
    IHSFinder(std::size_t snpLength, double cutoff, double minMAF, double scale, double binFactor);
    FreqVecMap unStdIHSByFreq() const { return m_unStandIHSByFreq; }
    LineMap    unStdIHSByLine() const { return m_unStandIHSByLine; }
    LineMap    freqsByLine() const    { return m_freqsByLine; }
    
    void run(HapMap* map, std::size_t start, std::size_t end);
    void runXpehh(HapMap* mA, HapMap* mB, std::size_t start, std::size_t end);
    LineMap normalize();
    
    void addData(const LineMap& freqsBySite, const LineMap& unStandIHSByLine, const FreqVecMap& unStandIHSbyLine);
    
protected:
    void processEHH(const EHH& ehh, std::size_t line);
    void processXPEHH(std::pair< EHH, EHH > e, size_t line);
    
    std::size_t m_snpLength;
    double m_cutoff;
    double m_minMAF;
    double m_scale;
    double m_binFactor;
    
    std::mutex m_mutex;
    std::mutex m_freqmutex;
    LineMap    m_freqsByLine;
    LineMap    m_unStandIHSByLine;
    FreqVecMap m_unStandIHSByFreq;
    LineMap    m_standIHSSingle;
    
    std::atomic<unsigned long long> m_counter;
};


#endif // IHSFINDER_H

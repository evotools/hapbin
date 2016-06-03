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

#include "ihsfinder.hpp"

IHSFinder::IHSFinder(std::size_t snpLength, double cutoff, double minMAF, double scale, unsigned long long maxExtend, int bins)
    : m_snpLength(snpLength), m_cutoff(cutoff), m_minMAF(minMAF), m_scale(scale), m_maxExtend(maxExtend), m_bins(bins), m_counter{}, m_reachedEnd{}, m_outsideMaf{}, m_nanResults{}
{}

void IHSFinder::processEHH(const EHH& ehh, std::size_t line)
{
    if (ehh.num + ehh.numNot != m_snpLength)
        return;

    double freqs = 0.0;
    double iHS = 0.0;
    if (ehh.iHH_0 > 0)
    {
        iHS = log(ehh.iHH_1/ehh.iHH_0);
        if (iHS == -std::numeric_limits<double>::infinity() || iHS == std::numeric_limits<double>::infinity())
        {
            ++m_nanResults;
            return;
        }
    } else {
        iHS = NAN;
        return;
    }

    freqs = ((int) (m_bins*ehh.num/(double)m_snpLength))/(double)m_bins;
    
    if (ehh.iHH_0 > 0)
    {
        m_freqmutex.lock();
        
        m_unStandIHSByFreq[freqs].push_back(iHS);
        m_freqmutex.unlock();
    }
    
    m_mutex.lock();
    m_freqsByLine[line] = freqs;
    m_unStandIHSByLine[line] = IhsScore(iHS, ehh.iHH_0, ehh.iHH_1);
    m_mutex.unlock();
}

void IHSFinder::processXPEHH(XPEHH& e, size_t line)
{
    if (e.iHH_A1 == 0.0 || e.iHH_B1 == 0.0)
        return;
    e.xpehh = log(e.iHH_A1/e.iHH_B1);
    m_mutex.lock();
    m_unStandXIHSByLine[line] = std::move(e);
    m_mutex.unlock();
}

IHSFinder::LineMap IHSFinder::normalize()
{
    StatsMap iHSStatsByFreq;
    
    for(const auto& it : m_unStandIHSByFreq)
    {
        iHSStatsByFreq[it.first] = stats(it.second);
    }
    
    for (const auto& it : m_unStandIHSByLine)
    {
        double freq = m_freqsByLine[it.first];
        if (iHSStatsByFreq[freq].stddev == 0)
            continue;
        m_standIHSSingle[it.first] = (it.second.iHS - iHSStatsByFreq[freq].mean)/iHSStatsByFreq[freq].stddev;
    }
    
    return m_standIHSSingle;
}

void IHSFinder::addData(const IHSFinder::LineMap& freqsBySite, 
                        const IHSFinder::IhsInfoMap& unStandIHSByLine, 
                        const IHSFinder::FreqVecMap& unStandIHSByFreq, 
                        unsigned long long reachedEnd, 
                        unsigned long long outsideMaf,
                        unsigned long long nanResults)
{
    m_reachedEnd += reachedEnd;
    m_outsideMaf += outsideMaf;
    m_nanResults += nanResults;
    m_mutex.lock();
    m_freqsByLine.insert(freqsBySite.begin(), freqsBySite.end());
    m_unStandIHSByLine.insert(unStandIHSByLine.begin(), unStandIHSByLine.end());
    m_mutex.unlock();
    m_freqmutex.lock();
    for (auto pair : unStandIHSByFreq)
    {
        std::vector<double>& v = m_unStandIHSByFreq[pair.first];
        v.insert(v.end(), pair.second.begin(), pair.second.end());
    }
    m_freqmutex.unlock();
}

void IHSFinder::addXData(const IHSFinder::LineMap& freqsBySite, 
                        const IHSFinder::XIhsInfoMap& unStandIHSByLine, 
                        const IHSFinder::FreqVecMap& unStandIHSByFreq, 
                        unsigned long long reachedEnd, 
                        unsigned long long outsideMaf,
                        unsigned long long nanResults)
{
    m_reachedEnd += reachedEnd;
    m_outsideMaf += outsideMaf;
    m_nanResults += nanResults;
    m_mutex.lock();
    m_freqsByLine.insert(freqsBySite.begin(), freqsBySite.end());
    m_unStandXIHSByLine.insert(unStandIHSByLine.begin(), unStandIHSByLine.end());
    m_mutex.unlock();
    m_freqmutex.lock();
    for (auto pair : unStandIHSByFreq)
    {
        std::vector<double>& v = m_unStandIHSByFreq[pair.first];
        v.insert(v.end(), pair.second.begin(), pair.second.end());
    }
    m_freqmutex.unlock();
}

void IHSFinder::runXpehh(HapMap* mA, HapMap* mB, std::size_t start, std::size_t end)
{
    #pragma omp parallel shared(mA,mB,start,end)
    {
        EHHFinder finder(mA->snpDataSize(), mB->snpDataSize(), 2000, m_cutoff, m_minMAF, m_scale, m_maxExtend);
        #pragma omp for schedule(dynamic,10)
        for(size_t i = start; i < end; ++i)
        {
            XPEHH xpehh = finder.findXPEHH(mA, mB, i, &m_reachedEnd);
            processXPEHH(xpehh, i);
            ++m_counter;
            unsigned long long tmp = m_counter;
            if (tmp % 1000 == 0)
            {
                std::cout << '\r' << tmp << "/" << (end-start);
            }
        }
    }
    std::cout << std::endl;
}

void IHSFinder::run(HapMap* map, std::size_t start, std::size_t end)
{
    #pragma omp parallel shared(map, start, end)
    {
        EHHFinder finder(map->snpDataSize(), map->snpDataSize(), 2000, m_cutoff, m_minMAF, m_scale, m_maxExtend);
        #pragma omp for schedule(dynamic,10)
        for(size_t i = start; i < end; ++i)
        {
            EHH ehh = finder.find(map, i, &m_reachedEnd, &m_outsideMaf);
            processEHH(ehh, i);
            ++m_counter;
            unsigned long long tmp = m_counter;
            if (tmp % 1000 == 0)
            {
                std::cout << '\r' << tmp << "/" << (end-start);
            }
        }
    }
    std::cout << std::endl;
}

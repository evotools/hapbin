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

#include "ihsfinder.hpp"

IHSFinder::IHSFinder(std::size_t snpLength,
                     double cutoff,
                     double minMAF,
                     double scale,
                     unsigned long long maxExtend,
                     int bins)
    : m_snpLength(snpLength)
    , m_cutoff(cutoff)
    , m_minMAF(minMAF)
    , m_scale(scale)
    , m_maxExtend(maxExtend)
    , m_bins(bins)
    , m_counter{}
    , m_reachedEnd{}
    , m_outsideMaf{}
    , m_nanResults{}

{}

void IHSFinder::processEHH(const EHH& ehh, std::size_t index)
{
    if ((ehh.num + ehh.numNot != m_snpLength) || (ehh.iHH_0 <= 0 && ehh.sl_0 <= 0))
        return;

    double freqs = ((int) (m_bins*ehh.num/(double)m_snpLength))/(double)m_bins;
    double ihs = 0.0, nsl = 0.0;
    if (ehh.iHH_0 > 0 || ehh.sl_0 > 0)
    {
        ihs = log(ehh.iHH_1/ehh.iHH_0);
        nsl = log(ehh.sl_1/ehh.sl_0);
        bool ihs_ok = (ihs != -std::numeric_limits<double>::infinity() && ihs != std::numeric_limits<double>::infinity() && ehh.iHH_0 > 0);
        bool nsl_ok = (nsl != -std::numeric_limits<double>::infinity() && nsl != std::numeric_limits<double>::infinity() && ehh.sl_0 > 0);
        if (ihs_ok || nsl_ok)
        {
            m_mutex.lock();
            if (ihs_ok)
                m_unStdIHSByFreq[freqs].push_back(ihs);
            if (nsl_ok)
                m_unStdNSLByFreq[freqs].push_back(nsl);
            m_unStdByIndex[index] = IhsScore(ihs, ehh.iHH_0, ehh.iHH_1, nsl, ehh.sl_0, ehh.sl_1);
            m_freqsByIndex[index] = freqs;
            m_mutex.unlock();
        }
        else
        {
            ++m_nanResults;
        }
    }
}

void IHSFinder::processXPEHH(XPEHH& e, size_t index)
{
    if (e.iHH_A1 == 0.0 || e.iHH_B1 == 0.0)
        return;
    e.xpehh = log(e.iHH_A1/e.iHH_B1);
    m_mutex.lock();
    m_unStdXIHSByIndex[index] = std::move(e);
    m_mutex.unlock();
}

void IHSFinder::normalize()
{
    StatsMap iHSStatsByFreq, nSLStatsByFreq;
    m_std.clear();

    for (const auto& it : m_unStdIHSByFreq)
        iHSStatsByFreq[it.first] = stats(it.second);
    for (const auto& it : m_unStdNSLByFreq)
        nSLStatsByFreq[it.first] = stats(it.second);

    for (const auto& it : m_unStdByIndex)
    {
        double freq = m_freqsByIndex[it.first];
        double stdIHS = NAN, stdNSL = NAN;
        if (iHSStatsByFreq[freq].stddev != 0)
            stdIHS = (it.second.iHS - iHSStatsByFreq[freq].mean)/iHSStatsByFreq[freq].stddev;
        if (nSLStatsByFreq[freq].stddev != 0)
            stdNSL = (it.second.nSL - nSLStatsByFreq[freq].mean)/nSLStatsByFreq[freq].stddev;
        IhsScore res = it.second;
        res.iHS = stdIHS;
        res.nSL = stdNSL;
        m_std[it.first] = res;
    }
}

void IHSFinder::addData(const IHSFinder::IndexMap& freqsByIndex,
                        const IHSFinder::IhsInfoMap& unStdIHSByIndex,
                        const IHSFinder::FreqVecMap& unStdIHSByFreq,
                        const IHSFinder::FreqVecMap& unStdNSLByFreq,
                        unsigned long long reachedEnd,
                        unsigned long long outsideMaf,
                        unsigned long long nanResults)
{
    m_reachedEnd += reachedEnd;
    m_outsideMaf += outsideMaf;
    m_nanResults += nanResults;
    m_mutex.lock();
    m_freqsByIndex.insert(freqsByIndex.begin(), freqsByIndex.end());
    m_unStdByIndex.insert(unStdIHSByIndex.begin(), unStdIHSByIndex.end());
    for (const auto& pair : unStdIHSByFreq)
    {
        auto& v = m_unStdIHSByFreq[pair.first];
        v.insert(v.end(), pair.second.cbegin(), pair.second.cend());
    }
    for (const auto& pair : unStdNSLByFreq)
    {
        auto& v = m_unStdNSLByFreq[pair.first];
        v.insert(v.end(), pair.second.cbegin(), pair.second.cend());
    }
    m_mutex.unlock();
}

void IHSFinder::addXData(const IHSFinder::IndexMap& freqsByIndex,
                        const IHSFinder::XIhsInfoMap& unStdXIHSByIndex,
                        const IHSFinder::FreqVecMap& unStdIHSByFreq,
                        unsigned long long reachedEnd,
                        unsigned long long outsideMaf,
                        unsigned long long nanResults)
{
    m_reachedEnd += reachedEnd;
    m_outsideMaf += outsideMaf;
    m_nanResults += nanResults;
    m_mutex.lock();
    m_freqsByIndex.insert(freqsByIndex.begin(), freqsByIndex.end());
    m_unStdXIHSByIndex.insert(unStdXIHSByIndex.begin(), unStdXIHSByIndex.end());
    for (auto pair : unStdIHSByFreq)
    {
        auto& v = m_unStdIHSByFreq[pair.first];
        v.insert(v.end(), pair.second.begin(), pair.second.end());
    }
    m_mutex.unlock();
}

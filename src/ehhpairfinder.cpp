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

#include "ehhpairfinder.hpp"
#include "ehh.hpp"
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cstdio>

EhhPairFinder::EhhPairFinder(const HapMap* hm, double cutoff, double minMAF, double scale, std::size_t maxBreadth, bool saveEhh)
    : m_hm(hm)
    , m_hd(hm->rawData())
    , m_parent{}
    , m_branch{}
    , m_parentNot{}
    , m_branchNot{}
    , m_ihh{}
    , m_ihhNot{}
    , m_sl{}
    , m_slNot{}
    , m_ehh{}
    , m_ehhNot{}
    , m_lastEhh{}
    , m_lastEhhNot{}
    , m_freq2{}
    , m_freqNot2{}
    , m_af{}
    , m_count{}
    , m_countNot{}
    , m_single{}
    , m_singleNot{}
    , m_calc{}
    , m_cutoff(cutoff)
    , m_minMAF(minMAF)
    , m_scale(scale)
    , m_maxBreadth(maxBreadth)
    , m_parentCount{}
    , m_branchCount{}
    , m_parentNotCount{}
    , m_branchNotCount{}
    , m_snpDataSize(hm->snpDataSize())
    , m_snpDataSizeULL(hm->snpDataSizeULL())
    , m_saveEhh(saveEhh)
#if VEC==4
    , m_mask(::bitsetMask4(hm->snpLength()))
#elif VEC==2
    , m_mask(::bitsetMask2(hm->snpLength()))
#else
    , m_mask(::bitsetMask<HapMap::PrimitiveType>(hmA->snpLength()))
#endif
    //, m_calced{}
{
    allocate();
}

void EhhPairFinder::allocate()
{
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        m_parent[i] = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSize*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branch[i] = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSize*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parentNot[i] = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSize*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branchNot[i] = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSize*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
}

void EhhPairFinder::reallocBranch()
{
    std::cout << "Reallocating branch buffer. new size: " << m_maxBreadth << " SNPs" << std::endl;
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        free(m_branch[i]);
        free(m_branchNot[i]);
        m_branch[i] = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSize*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_branchNot[i] = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSize*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
}

void EhhPairFinder::setCore(std::size_t focus1, std::size_t focus2)
{
    for (std::size_t i = 0; i < m_snpDataSize; ++i)
    {
        HapMap::PrimitiveType b00 = ~m_hd[focus1*m_snpDataSize+i] & ~m_hd[focus2*m_snpDataSize+i];
        HapMap::PrimitiveType b01 = ~m_hd[focus1*m_snpDataSize+i] &  m_hd[focus2*m_snpDataSize+i];
        HapMap::PrimitiveType b10 =  m_hd[focus1*m_snpDataSize+i] & ~m_hd[focus2*m_snpDataSize+i];
        HapMap::PrimitiveType b11 =  m_hd[focus1*m_snpDataSize+i] &  m_hd[focus2*m_snpDataSize+i];

        m_parent[0b00][i] = b00;
        m_parent[0b01][i] = b01;
        m_parent[0b10][i] = b10;
        m_parent[0b11][i] = b11;

        m_parentNot[0b00][i] = ~b00;
        m_parentNot[0b01][i] = ~b01;
        m_parentNot[0b10][i] = ~b10;
        m_parentNot[0b11][i] = ~b11;
    }
    m_parent   [0b00][m_snpDataSize-1] &= m_mask;
    m_parentNot[0b01][m_snpDataSize-1] &= m_mask;
    m_parentNot[0b10][m_snpDataSize-1] &= m_mask;
    m_parentNot[0b11][m_snpDataSize-1] &= m_mask;
    m_parentCount[0b00] = 1UL;
    m_parentCount[0b01] = 1UL;
    m_parentCount[0b10] = 1UL;
    m_parentCount[0b11] = 1UL;
    m_parentNotCount[0b00] = 1UL;
    m_parentNotCount[0b01] = 1UL;
    m_parentNotCount[0b10] = 1UL;
    m_parentNotCount[0b11] = 1UL;
}

int EhhPairFinder::countHaps(const HapMap::PrimitiveType *haps)
{
    int count = 0;
    unsigned long long* hapsULL = (unsigned long long*) haps;
    for(std::size_t i = 0; i < m_snpDataSizeULL; ++i)
    {
        count += popcount1(hapsULL[i]);
    }
    return count;
}

bool EhhPairFinder::branch(std::size_t core, std::size_t index, std::size_t bufferOffset, HapMap::PrimitiveType *parent[4], HapMap::PrimitiveType *branch[4], std::size_t (&branchCount)[4], double freq2[4], double (&ehh)[4], int (&single)[4])
{
    int count = countHaps(&parent[core][bufferOffset*m_snpDataSize]);
    if (count == 0)
    {
        return false;
    }
    else if (count == 1)
    {
        ++single[core];
    }
    else
    {
        ehh[core] += ((double)count*count)*(freq2[core]);

        for(std::size_t j = 0; j < m_snpDataSize; ++j)
        {
            branch[core][branchCount[core]*m_snpDataSize+j] = parent[core][bufferOffset*m_snpDataSize+j] & m_hd[index*m_snpDataSize+j];
        }
        ++branchCount[core];
        for(std::size_t j = 0; j < m_snpDataSize-1; ++j)
        {
            branch[core][branchCount[core]*m_snpDataSize+j] = parent[core][bufferOffset*m_snpDataSize+j] & ~m_hd[index*m_snpDataSize+j];
        }
        branch[core][branchCount[core]*m_snpDataSize+m_snpDataSize-1] = parent[core][bufferOffset*m_snpDataSize+m_snpDataSize-1] & ~m_hd[index*m_snpDataSize+m_snpDataSize-1] & m_mask;
        ++branchCount[core];
    }
    if (branchCount[core] > m_maxBreadth - 2)
    {
        return true;
    }
    return false;
}

void EhhPairFinder::createBranches(std::size_t index)
{
    bool overflow;
    bool realloced = false;
    do
    {
        overflow = false;
        for(std::size_t i = 0b00; i <= 0b11; ++i)
        {
            overflow = processBuffer(i, index);
            if (overflow)
                break;
        }
        if (overflow)
        {
            realloced = true;
            m_maxBreadth += 100;
            reallocBranch();
        }
    }
    while (overflow);
    if (realloced)
        reallocBranch();

}

bool EhhPairFinder::processBuffer(std::size_t core, std::size_t index)
{
    bool overflow;
    m_ehh[core] = 0.0;
    m_ehhNot[core] = 0.0;
    m_branchCount[core] = 0UL;
    m_branchNotCount[core] = 0UL;
    int single[4]{};
    int singleNot[4]{};
    if (!m_calc[core])
        return false;
    for(std::size_t i = 0; i < m_parentCount[core]; ++i)
    {
        overflow = branch(core, index, i, m_parent, m_branch, m_branchCount, m_freq2, m_ehh, single);
        if (overflow)
            return true;
    }
    for(std::size_t i = 0; i < m_parentNotCount[core]; ++i)
    {
        overflow = branch(core, index, i, m_parentNot, m_branchNot, m_branchNotCount, m_freqNot2, m_ehhNot, singleNot);
        if (overflow)
            return true;
    }
    std::swap(m_parent[core], m_branch[core]);
    std::swap(m_parentNot[core], m_branchNot[core]);
    for (std::size_t core = 0b00; core <= 0b11; ++core)
    {
        m_single[core] += single[core];
    }
    for (std::size_t core = 0b00; core <= 0b11; ++core)
    {
        m_singleNot[core] += singleNot[core];
    }

    //std::cout << core << " single: " << m_single[core] << " " << ((double)m_single[core])*m_freq2[core] << std::endl;
    m_ehh[core] += ((double)m_single[core])*m_freq2[core];
    m_ehhNot[core] += ((double)m_singleNot[core])*m_freqNot2[core];

    m_parentCount[core] = m_branchCount[core];
    m_parentNotCount[core] = m_branchNotCount[core];

    return false;
}

void EhhPairFinder::calcBranch(std::size_t core, std::size_t index, bool right)
{
    std::size_t offset1, offset2;
    if (!right)
    {
        offset1 = 2;
        offset2 = 1;
    }
    else
    {
        offset1 = -1;
        offset2 = -2;
    }
    double scale = (double)(m_scale) / (double)(m_hm->physicalPosition(index+offset1) - m_hm->physicalPosition(index+offset2));
    if (scale > 1)
        scale=1;
    m_ihh[core] += (m_hm->geneticPosition(index+offset1)-m_hm->geneticPosition(index+offset2))*(m_lastEhh[core] + m_ehh[core])*scale*0.5;
    m_ihhNot[core] += (m_hm->geneticPosition(index+offset1)-m_hm->geneticPosition(index+offset2))*(m_lastEhhNot[core] + m_ehhNot[core])*scale*0.5;

    if (m_ehh[core] < m_lastEhh[core])
        m_sl[core] += m_ehh[core];
    if (m_ehhNot[core] < m_lastEhhNot[core])
        m_slNot[core] += m_ehh[core];

    m_lastEhh[core] = m_ehh[core];
    m_lastEhhNot[core] = m_ehhNot[core];
}

bool EhhPairFinder::calcBranches(std::size_t index, bool right)
{
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        if ((m_lastEhh[i] <= m_cutoff + 1e-15 && m_lastEhhNot[i] <= m_cutoff + 1e-15) || !m_calc[i])
            continue;
        calcBranch(i, index, right);
    }
}

void EhhPairFinder::calcFreq()
{
    for(size_t i = 0b00; i <= 0b11; ++i)
    {
        int count = countHaps(m_parent[i]);
        int countNot = m_hm->snpLength() - count;
        int mac = (count > countNot) ? countNot : count;
        m_count[i] = count;
        m_countNot[i] = countNot;
        m_af[i] = (double)count/(double)m_hm->snpLength();
        double maf = (double)mac/(double)m_hm->snpLength();
        m_calc[i] = (maf > m_minMAF);
        //std::cout << m_calc[i] << " ";
        m_freq2[i] = 1.0/(double)(count*count);
        m_freqNot2[i] = 1.0/(double)(countNot*countNot);
    }
    //std::cout << std::endl;
}

bool EhhPairFinder::passesAf(std::size_t focus1, std::size_t focus2)
{
    setCore(focus1, focus2);
    calcFreq();

    int numCalc = 0;
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        if (m_calc[i])
            numCalc++;
    }
    if (numCalc >= 3)
        return true;
    else
        return false;
}

bool EhhPairFinder::reachedCutoff(double ehh)
{
    return ((ehh <= m_cutoff + 1e-15) || std::isinf(ehh) || std::isnan(ehh));
}

bool EhhPairFinder::reachedCutoffs()
{
    return ((reachedCutoff(m_lastEhh[0b00]) || reachedUniqueness(0b00)) &&
            (reachedCutoff(m_lastEhh[0b01]) || reachedUniqueness(0b01)) &&
            (reachedCutoff(m_lastEhh[0b10]) || reachedUniqueness(0b10)) &&
            (reachedCutoff(m_lastEhh[0b11]) || reachedUniqueness(0b11)) &&
            (reachedCutoff(m_lastEhhNot[0b00]) || reachedUniqueness(0b00)) &&
            (reachedCutoff(m_lastEhhNot[0b01]) || reachedUniqueness(0b01)) &&
            (reachedCutoff(m_lastEhhNot[0b10]) || reachedUniqueness(0b10)) &&
            (reachedCutoff(m_lastEhhNot[0b11]) || reachedUniqueness(0b11)));
}

bool EhhPairFinder::reachedUniqueness(std::size_t core)
{
    return ((m_single[core] == m_count[core]) || (m_singleNot[core] == m_countNot[core]));
}

bool EhhPairFinder::reachedUniqueness()
{
    return (reachedUniqueness(0b00) &&
            reachedUniqueness(0b01) &&
            reachedUniqueness(0b10) &&
            reachedUniqueness(0b11));
}

bool EhhPairFinder::checkMaf(std::size_t focus1, std::size_t focus2)
{
    int count1 = countHaps(&m_hd[focus1*m_snpDataSize]);
    int count2 = countHaps(&m_hd[focus2*m_snpDataSize]);
    int count1Not = m_hm->snpLength() - count1;
    int count2Not = m_hm->snpLength() - count2;
    int count1Min = (count1 > count1Not)? count1Not : count1;
    int count2Min = (count2 > count2Not)? count2Not : count2;
    double maf1 = (double)count1Min/(double)m_hm->snpLength();
    double maf2 = (double)count2Min/(double)m_hm->snpLength();
    //std::cout << "maf: " << maf1 << " " << maf2 << std::endl;
    return (maf1 >= m_minMAF && maf2 >= m_minMAF);
}

EHHPair EhhPairFinder::calcEhhPair(std::size_t focus1, std::size_t focus2, bool focus1Conditioned)
{
    EHHPair ret{};
    if (!checkMaf(focus1, focus2))
        return ret;
    if (focus1 < 2)
        return ret;

    setCore(focus1, focus2);
    calcFreq();
    createBranches(focus1-1);
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        m_lastEhh[i] = (m_calc[i])? 1.0 : 0.0;
        m_lastEhhNot[i] = (m_calc[i])? 1.0 : 0.0;
        m_ihh[i] = 0.0;
        m_ihhNot[i] = 0.0;
        m_sl[i] = 0.0;
        m_slNot[i] = 0.0;
        m_single[i] = 0;
        m_singleNot[i] = 0;
    }
    for (std::size_t currIndex = focus1 - 2;; --currIndex)
    {
        //std::cout << "A " << m_lastEhh[0] << " " << m_lastEhh[1] << " " << m_lastEhh[2] << " " << m_lastEhh[3] << " " << m_lastEhhNot[0] << " " << m_lastEhhNot[1] << " " << m_lastEhhNot[2] << " " << m_lastEhhNot[3] << std::endl;//" " << m_ihh[0] << " " << m_ihh[3] << " " << m_ihhNot[0] << " " << m_ihhNot[3]  << std::endl;
        if (m_saveEhh)
            ret.upstream.push_back(HapPairStats(m_lastEhh, m_lastEhhNot));
        if (reachedCutoffs())
            break;
        createBranches(currIndex);
        calcBranches(currIndex, false);
        if (reachedUniqueness())
            break;
        if (currIndex == 0)
        {
            //std::cout << "reached start of chromosome. " << m_lastEhh[0] << " " << m_lastEhh[1] << " " << m_lastEhh[2] << " " << m_lastEhh[3] << " "
            //                                             << m_lastEhhNot[0] << " " << m_lastEhhNot[1] << " " << m_lastEhhNot[2] << " " << m_lastEhh[3] << " " << std::endl;
            return ret;
        }
    }

    setCore(focus1, focus2);
    std::size_t rightFocus;
    if (focus1Conditioned)
        rightFocus = focus1;
    else
        rightFocus = focus2;
    createBranches(rightFocus+1);
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        m_lastEhh[i] = (m_calc[i])? 1.0 : 0.0;
        m_lastEhhNot[i] = (m_calc[i])? 1.0 : 0.0;
        m_single[i] = 0;
        m_singleNot[i] = 0;
    }
    for (std::size_t currIndex = rightFocus + 2; currIndex < m_hm->numSnps(); ++currIndex)
    {
        //std::cout << "B " << m_lastEhh[0] << " " << m_lastEhh[1] << " " << m_lastEhh[2] << " " << m_lastEhh[3] << " " << m_lastEhhNot[0] << " " << m_lastEhhNot[1] << " " << m_lastEhhNot[2] << " " << m_lastEhhNot[3] << std::endl;//" " << m_ihh[0] << " " << m_ihh[3] << " " << m_ihhNot[0] << " " << m_ihhNot[3]  << std::endl;
        if (m_saveEhh)
            ret.downstream.push_back(HapPairStats(m_lastEhh, m_lastEhhNot));
        if (reachedCutoffs())
            break;
        createBranches(currIndex);
        calcBranches(currIndex, true);
        if (reachedUniqueness())
            break;
        if (currIndex == m_hm->numSnps()-1)
        {
            //std::cout << "reached end of chromosome." << std::endl;
            return ret;
        }
    }
    ret.focus[0] = focus1;
    ret.focus[1] = focus2;
    std::copy_n(m_count, 4, ret.count);
    std::copy_n(m_af, 4, ret.af);
    std::copy_n(m_ihh, 4, ret.ihh);
    std::copy_n(m_ihhNot, 4, ret.ihhNot);
    std::copy_n(m_sl, 4, ret.sl);
    std::copy_n(m_slNot, 4, ret.slNot);
    //std::cout << "iHH:    " << m_ihh[0b00] << " " << m_ihh[0b01] << " " << m_ihh[0b10] << " " << m_ihh[0b11] << std::endl;
    //std::cout << "iHHNot: " << m_ihhNot[0b00] << " " << m_ihhNot[0b01] << " " << m_ihhNot[0b10] << " " << m_ihhNot[0b11] << std::endl;
    //std::cout << "iHS: " << log(m_ihh[0b00]/m_ihhNot[0b00]) << " " << log(m_ihh[0b01]/m_ihhNot[0b01]) << " " << log(m_ihh[0b10]/m_ihhNot[0b10]) << " " << log(m_ihh[0b11]/m_ihhNot[0b11]) << std::endl;
    return ret;
}

EhhPairFinder::~EhhPairFinder()
{
    for(std::size_t i = 0b00; i <= 0b11; ++i)
    {
        aligned_free(m_parent[i]);
        aligned_free(m_branch[i]);
        aligned_free(m_parentNot[i]);
        aligned_free(m_branchNot[i]);
    }
}

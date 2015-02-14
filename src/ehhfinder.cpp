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

#include "ehhfinder.hpp"
#include "hapmap.hpp"
#include <algorithm>

EHHFinder::EHHFinder(std::size_t snpDataSizeA, std::size_t snpDataSizeB, std::size_t maxBreadth, double cutoff, double minMAF, double scale)
    : m_maxSnpDataSize(snpDataSizeA)
    , m_maxBreadth(maxBreadth)
    , m_bufferSize(snpDataSizeA*maxBreadth)
    , m_cutoff(cutoff)
    , m_minMAF(minMAF)
    , m_scale(scale)
    , m_freqA{}
    , m_freqB{}
    , m_freqP{}
    , m_ehhA{}
    , m_ehhB{}
    , m_ehhP{}
{
    m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
}

/**
 * Set initial state. Set m_parent0 to '0' core haplotype positions, m_parent1 to '1' core haplotype positions.
 * 
 * calcBranch counts the state of the parent (previous) branch, not the counts of the new (current) level. Therefore, we must advance
 * the calculations by one iteration before starting. 
 */
void EHHFinder::setInitial(std::size_t focus, std::size_t line)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    m_single0count = 0ULL;
    m_single1count = 0ULL;
    

    for (int j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[               j] = ~m_hdA[focus*m_snpDataSizeA+j] &  m_hdA[line*m_snpDataSizeA+j];
        m_parent0[m_snpDataSizeA+j] = (~m_hdA[focus*m_snpDataSizeA+j]) & (~m_hdA[line*m_snpDataSizeA+j]);
        m_parent1[               j] =  m_hdA[focus*m_snpDataSizeA+j] &  m_hdA[line*m_snpDataSizeA+j];
        m_parent1[m_snpDataSizeA+j] =  m_hdA[focus*m_snpDataSizeA+j] & ~m_hdA[line*m_snpDataSizeA+j];
    }
    m_parent0[  m_snpDataSizeA-1] &= m_maskA;
    m_parent0[2*m_snpDataSizeA-1] &= m_maskA;
    m_parent1[2*m_snpDataSizeA-1] &= m_maskA;
}

void EHHFinder::setInitialXPEHH(std::size_t focus)
{
    for (std::size_t j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[j] = ~m_hdA[focus*m_snpDataSizeA+j];
    }
    m_parent0[m_snpDataSizeA-1] &= m_maskA;
    for (std::size_t j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[j+m_snpDataSizeA] = m_hdA[focus*m_snpDataSizeA+j];
    }
    for (std::size_t j = 0; j < m_snpDataSizeB; ++j)
    {
        m_parent1[j] = ~m_hdB[focus*m_snpDataSizeB+j];
    }
    m_parent1[m_snpDataSizeB-1] &= m_maskB;
    for (std::size_t j = 0; j < m_snpDataSizeB; ++j)
    {
        m_parent1[j+m_snpDataSizeB] = m_hdB[focus*m_snpDataSizeB+j];
    }
}

void EHHFinder::calcBranch(HapMap* hm, std::size_t focus, HapMap::PrimitiveType* parent, std::size_t parentcount, HapMap::PrimitiveType* branch, std::size_t& branchcount, std::size_t currLine, double freq, double &probs, std::size_t& singlecount, bool* overflow)
{
    HapMap::PrimitiveType* mapData = hm->rawData();
    std::size_t snpDataSize = hm->snpDataSize();
    std::size_t snpDataSizeULL = hm->snpDataSizeULL();
    std::size_t bcnt = 0;
    for(std::size_t i = 0; i < parentcount; ++i)
    {
        int count = 0;
        unsigned long long *leaf = (unsigned long long*) &parent[i*snpDataSize];
        for (std::size_t j = 0; j < snpDataSizeULL; ++j)
        {
            count += popcount1(leaf[j]);
        }
        
        if (count == 0)
        {
            continue;
        }
        else if (count == 1)
        {
            ++singlecount;
            continue;
        }
        else
        {
            probs += (count*freq)*(count*freq);
            for(std::size_t j = 0; j < snpDataSize; ++j)
            {
                branch[bcnt*snpDataSize+j] = parent[i*snpDataSize+j] & mapData[currLine*snpDataSize+j];
            }
            ++bcnt;
            for(std::size_t j = 0; j < snpDataSize-1; ++j)
            {
                branch[bcnt*snpDataSize+j] = parent[i*snpDataSize+j] & ~mapData[currLine*snpDataSize+j];
            }
            branch[bcnt*snpDataSize+snpDataSize-1] = (parent[i*snpDataSize+snpDataSize-1] & ~mapData[currLine*snpDataSize+snpDataSize-1]) & m_maskA;
            ++bcnt;
        }
        if (bcnt > m_maxBreadth-2)
        {
            *overflow = true;
            return;
        }
    }
    branchcount = bcnt;
}

void EHHFinder::calcBranchXPEHH(HapMap* hmA, HapMap* hmB, std::size_t currLine, std::size_t& single0,  std::size_t& single1, bool* overflow)
{
    for(std::size_t i = 0; i < m_parent0count; ++i)
    {
        int numA = 0;
        unsigned long long *parentA_ULL = (unsigned long long*) &m_parent0[i*m_snpDataSizeA];
        for (std::size_t j = 0; j < m_snpDataSizeULL_A; ++j)
        {
            numA += popcount1(parentA_ULL[j]);
        }
        int numB = 0;
        unsigned long long *parentB_ULL = (unsigned long long*) &m_parent1[i*m_snpDataSizeB];
        for (size_t j = 0; j < m_snpDataSizeULL_B; ++j)
        {
            numB += popcount1(parentB_ULL[j]);
        }
        int numPooled = numA + numB;
        //std::cout << i << "numA: " << numA << " numB: " << numB << " numP: " << numPooled << " freqP: " << m_freqP << " " << m_snpDataSizeUL_A  << std::endl;
        if (numPooled == 0)
        {
            continue;
        } 
        else if (numPooled == 1)
        {
            single0 += numA;
            single1 += numB;
            continue;
        } 
        else
        {
            m_ehhA += (numA*m_freqA)*(numA*m_freqA);
            m_ehhB += (numB*m_freqB)*(numB*m_freqB);
            m_ehhP += (numPooled*m_freqP)*(numPooled*m_freqP);
            
            //Population A
            for(std::size_t j = 0; j < m_snpDataSizeA; ++j)
            {
                m_branch0[m_branch0count*m_snpDataSizeA+j] = m_parent0[i*m_snpDataSizeA+j] & m_hdA[currLine*m_snpDataSizeA+j];
            }
            ++m_branch0count;
            for(std::size_t j = 0; j < m_snpDataSizeA-1; ++j)
            {
                m_branch0[m_branch0count*m_snpDataSizeA+j] = m_parent0[i*m_snpDataSizeA+j] & ~m_hdA[currLine*m_snpDataSizeA+j];
            }
            m_branch0[m_branch0count*m_snpDataSizeA+m_snpDataSizeA-1] = (m_parent0[i*m_snpDataSizeA+m_snpDataSizeA-1] & ~m_hdA[currLine*m_snpDataSizeA+m_snpDataSizeA-1]) & m_maskA;
            ++m_branch0count;
            
            //Population B
            for(std::size_t j = 0; j < m_snpDataSizeB; ++j)
            {
                m_branch1[m_branch1count*m_snpDataSizeB+j] = m_parent1[i*m_snpDataSizeB+j] & m_hdB[currLine*m_snpDataSizeB+j];
            }
            ++m_branch1count;
            for(std::size_t j = 0; j < m_snpDataSizeB-1; ++j)
            {
                m_branch1[m_branch1count*m_snpDataSizeB+j] = m_parent1[i*m_snpDataSizeB+j] & ~m_hdB[currLine*m_snpDataSizeB+j];
            }
            m_branch1[m_branch1count*m_snpDataSizeB+m_snpDataSizeB-1] = (m_parent1[i*m_snpDataSizeB+m_snpDataSizeB-1] & ~m_hdB[currLine*m_snpDataSizeB+m_snpDataSizeB-1]) & m_maskB;
            ++m_branch1count;
        }
        if (m_branch0count > m_maxBreadth-2)
        {
            *overflow = true;
            return;
        }
    }
}

void EHHFinder::calcBranchesXPEHH(std::size_t currLine)
{
    bool overflow;
    bool realloced = false;
    std::size_t single0{}, single1{};
    m_ehhA = 0.0;
    m_ehhB = 0.0;
    m_ehhP = 0.0;
    while (true)
    {
        overflow = false;
        single0 = 0; single1 = 0;
        m_branch0count = 0;
        m_branch1count = 0;
        calcBranchXPEHH(m_hmA, m_hmB, currLine, single0, single1, &overflow);
        if (overflow)
        {
            m_maxBreadth += 500;
            aligned_free(m_branch0);
            aligned_free(m_branch1);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            realloced = true;
        }
        else
            break;
    }
    if (realloced)
    {
        aligned_free(m_parent0);
        aligned_free(m_parent1);
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
    m_single0count += single0;
    m_single1count += single1;
    m_parent0count = m_branch0count;
    m_parent1count = m_branch1count;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    std::swap(m_parent0, m_branch0);
    std::swap(m_parent1, m_branch1);
}

std::pair< EHH, EHH > EHHFinder::findXPEHH(HapMap* hmA, HapMap* hmB, std::size_t focus)
{
    m_hmA = hmA;
    m_hmB = hmB;
    m_hdA = hmA->rawData();
    m_snpDataSizeA = hmA->snpDataSize();
    m_snpDataSizeULL_A = hmA->snpDataSizeULL();
    m_hdB = hmB->rawData();
    m_snpDataSizeB = hmB->snpDataSize();
    m_snpDataSizeULL_B = hmB->snpDataSizeULL();
#if VEC==4
    m_maskA = ::bitsetMask4(hmA->snpLength());
    m_maskB = ::bitsetMask4(hmB->snpLength());
#elif VEC==2
    m_maskA = ::bitsetMask2(hmA->snpLength());
    m_maskB = ::bitsetMask2(hmB->snpLength());
#else
    m_maskA = ::bitsetMask<HapMap::PrimitiveType>(hmA->snpLength());
    m_maskB = ::bitsetMask<HapMap::PrimitiveType>(hmB->snpLength());
#endif
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    m_single0count = 0ULL;
    m_single1count = 0ULL;
    std::pair<EHH,EHH> ret;
    ret.first.index = focus;
    for(std::size_t i = 0; i < m_snpDataSizeA; ++i)
    {
        ret.first.num += POPCOUNT(m_hdA[focus*m_snpDataSizeA+i]);
    }
    ret.first.numNot = hmA->snpLength() - ret.first.num;
    for(std::size_t i = 0; i < m_snpDataSizeB; ++i)
    {
        ret.second.num += POPCOUNT(m_hdB[focus*m_snpDataSizeB+i]);
    }
    ret.second.numNot = hmA->snpLength() - ret.first.num;
    m_freqA = 1.0/(double)hmA->snpLength();
    m_freqB = 1.0/(double)hmB->snpLength();
    m_freqP = 1.0/(double)(hmA->snpLength()+hmB->snpLength());
    double probASingle = m_freqA*m_freqA;
    double probBSingle = m_freqB*m_freqB;
    double probPSingle = m_freqP*m_freqP;
    double f = ret.first.num*m_freqA;
    double lastEhhA = f*f+(1.0-f)*(1.0-f);
    f = ret.second.num*m_freqB;
    double lastEhhB = f*f+(1.0-f)*(1.0-f);
    f = (ret.first.num+ret.second.num)*m_freqP;
    double lastEhhP = f*f+(1.0-f)*(1.0-f);
    
    setInitialXPEHH(focus);
    calcBranchesXPEHH(focus-1);
    m_single0count = 0ULL;
    m_single1count = 0ULL;
    if (focus > 1)
    {
        for (std::size_t currLine = focus - 2;; --currLine)
        {
            
            double scale = (double)(m_scale) / (double)(hmA->physicalPosition(currLine+2) - hmA->physicalPosition(currLine+1));
            if (scale > 1)
                scale=1;
            
            calcBranchesXPEHH(currLine);
            
            m_ehhA += probASingle*m_single0count;
            m_ehhB += probBSingle*m_single1count;
            m_ehhP += probPSingle*(m_single0count+m_single1count);
            
            if (lastEhhP <= m_cutoff + 1e-15)
                break;
            ret.first.iHH_d  += (hmA->geneticPosition(currLine+2)-hmA->geneticPosition(currLine+1))*(lastEhhA + m_ehhA)*scale*0.5;    
            ret.second.iHH_d += (hmA->geneticPosition(currLine+2)-hmA->geneticPosition(currLine+1))*(lastEhhB + m_ehhB)*scale*0.5;
            
            lastEhhA = m_ehhA;
            lastEhhB = m_ehhB;
            lastEhhP = m_ehhP;
            
            if ((m_single0count+m_single1count) == (hmA->snpLength()+hmB->snpLength()))
                break;
            if(currLine == 0)
                return std::pair<EHH,EHH>();
        }
    }
    
    f = ret.first.num*m_freqA;
    lastEhhA = f*f+(1.0-f)*(1.0-f);
    f = ret.second.num*m_freqB;
    lastEhhB = f*f+(1.0-f)*(1.0-f);
    f = (ret.first.num+ret.second.num)*m_freqP;
    lastEhhP = f*f+(1.0-f)*(1.0-f);
    
    setInitialXPEHH(focus);
    calcBranchesXPEHH(focus+1);
    m_single0count = 0ULL;
    m_single1count = 0ULL;
    for (std::size_t currLine = focus + 2; currLine < hmA->numSnps(); ++currLine)
    {
        double scale = (double)(m_scale) / (double)(hmA->physicalPosition(currLine-1) - hmA->physicalPosition(currLine-2));
        if (scale > 1)
            scale=1;
        
        calcBranchesXPEHH(currLine);
        
        m_ehhA += probASingle*m_single0count;
        m_ehhB += probBSingle*m_single1count;
        m_ehhP += probPSingle*(m_single0count+m_single1count);
        
        if (lastEhhP <= m_cutoff + 1e-15)
            break;
        ret.first.iHH_d  += (hmA->geneticPosition(currLine-1)-hmA->geneticPosition(currLine-2))*(lastEhhA + m_ehhA)*scale*0.5;    
        ret.second.iHH_d += (hmA->geneticPosition(currLine-1)-hmA->geneticPosition(currLine-2))*(lastEhhB + m_ehhB)*scale*0.5;
        
        lastEhhA = m_ehhA;
        lastEhhB = m_ehhB;
        lastEhhP = m_ehhP;
        
        if ((m_single0count+m_single1count) == (hmA->snpLength()+hmB->snpLength()))
            break;
        if (currLine == hmA->numSnps()-1)
            return std::pair<EHH,EHH>();
    }
    return ret;
}

void EHHFinder::calcBranches(HapMap* hapmap, std::size_t focus, std::size_t currLine, double freq0,  double freq1, HapStats &stats)
{
    bool overflow;
    bool realloced = false;
    std::size_t single0{};
    std::size_t single1{};
    while(true)
    {
        overflow = false;
        single0 = 0; single1 = 0;
        calcBranch(hapmap, focus, m_parent0, m_parent0count, m_branch0, m_branch0count, currLine, freq0, stats.probsNot, single0, &overflow);
        if(overflow)
        {
            m_maxBreadth += 100;
            aligned_free(m_branch0);
            aligned_free(m_branch1);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            realloced = true;
            continue;
        }
        overflow = false;
        calcBranch(hapmap, focus, m_parent1, m_parent1count, m_branch1, m_branch1count, currLine, freq1, stats.probs, single1, &overflow);
        if(overflow)
        {
            m_maxBreadth += 100;
            aligned_free(m_branch0);
            aligned_free(m_branch1);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
            realloced = true;
            continue;
        }
        break;
    }
    if (realloced)
    {
        aligned_free(m_parent0);
        aligned_free(m_parent1);
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth*sizeof(HapMap::PrimitiveType)));
    }
    m_parent0count = m_branch0count;
    m_parent1count = m_branch1count;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    m_single0count += single0;
    m_single1count += single1;
    std::swap(m_parent0, m_branch0);
    std::swap(m_parent1, m_branch1);
}

EHH EHHFinder::find(HapMap* hapmap, std::size_t focus, bool ehhsave)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    m_single0count = 0ULL;
    m_single1count = 0ULL;
    m_hdA = hapmap->rawData();
    m_snpDataSizeA = hapmap->snpDataSize();
#if VEC==4
    m_maskA = ::bitsetMask4(hapmap->snpLength());
#elif VEC==2
    m_maskA = ::bitsetMask2(hapmap->snpLength());
#else
    m_maskA = ::bitsetMask<HapMap::PrimitiveType>(hapmap->snpLength());
#endif
    EHH ret;
    ret.index = focus;
    
    for(std::size_t i = 0; i < m_snpDataSizeA; ++i)
    {
        ret.num += POPCOUNT(m_hdA[focus*m_snpDataSizeA+i]);
    }
    ret.numNot = hapmap->snpLength() - ret.num;
    
    double maxEHH = ret.num/(double)hapmap->snpLength();
    if (!(maxEHH <= 1.0 - m_minMAF && maxEHH >= m_minMAF)) {
        return EHH();
    }
    if (focus < 2)
        return EHH();
    
    double freq0 = 1.0/(double)ret.numNot;
    double freq1 = 1.0/(double)ret.num;
    double probSingle = freq1*freq1;
    double probNotSingle = freq0*freq0;
    double probs, probsNot;
    double lastProbs = 1.0, lastProbsNot = 1.0;
    
    setInitial(focus, focus-1);
    
    for (std::size_t currLine = focus - 2;; --currLine)
    {
        HapStats stats;
        double scale = (double)(m_scale) / (double)(hapmap->physicalPosition(currLine+2) - hapmap->physicalPosition(currLine+1));
        if (scale > 1)
            scale=1;
        
        calcBranches(hapmap, focus, currLine, freq0, freq1, stats);
        
        stats.probs +=  probSingle*m_single1count;
        stats.probsNot += probNotSingle*m_single0count;

        if (lastProbs > m_cutoff + 1e-15)
            ret.iHH_d += (hapmap->geneticPosition(currLine+2)-hapmap->geneticPosition(currLine+1))*(lastProbs + stats.probs)*scale*0.5;    
        if (lastProbsNot > m_cutoff + 1e-15)
            ret.iHH_a += (hapmap->geneticPosition(currLine+2)-hapmap->geneticPosition(currLine+1))*(lastProbsNot + stats.probsNot)*scale*0.5;   
        
        lastProbs = stats.probs;
        lastProbsNot = stats.probsNot;
        if (ehhsave)
            ret.upstream.push_back(std::move(stats));
        
        
        if (lastProbs <= m_cutoff + 1e-15 && lastProbsNot <= m_cutoff + 1e-15)
            break;
        if ((m_single0count+m_single1count) == hapmap->snpLength())
            break;
        if (currLine == 0)
        {
            return EHH();
        }
    }
    
    setInitial(focus,focus+1);
    lastProbs = 1.0, lastProbsNot = 1.0;
    for (std::size_t currLine = focus + 2; currLine < hapmap->numSnps(); ++currLine)
    {
        HapStats stats;
        double scale = double(m_scale) / double(hapmap->physicalPosition(currLine-1) - hapmap->physicalPosition(currLine-2));
        if (scale > 1)
            scale=1;
        
        int core0 = 0, core1 = 0;
        
        calcBranches(hapmap, focus, currLine, freq0, freq1, stats);
        
        stats.probs +=  probSingle*m_single1count;
        stats.probsNot += probNotSingle*m_single0count;
        
        if (lastProbs > m_cutoff + 1e-15) {
            ret.iHH_d += (hapmap->geneticPosition(currLine-1)-hapmap->geneticPosition(currLine-2))*(lastProbs + stats.probs)*scale*0.5;
        }
        if (lastProbsNot > m_cutoff + 1e-15) {
            ret.iHH_a += (hapmap->geneticPosition(currLine-1)-hapmap->geneticPosition(currLine-2))*(lastProbsNot + stats.probsNot)*scale*0.5;
        }
        
        lastProbs = stats.probs;
        lastProbsNot = stats.probsNot;
        if (ehhsave)
            ret.downstream.push_back(std::move(stats));
        
        if (lastProbs <= m_cutoff + 1e-15 && lastProbsNot <= m_cutoff + 1e-15 || (m_single0count+m_single1count) == hapmap->snpLength())
            break;
        
        if (currLine == hapmap->numSnps()-1)
            return EHH();
    }
    return ret;
}

EHHFinder::~EHHFinder()
{
    aligned_free(m_branch0);
    aligned_free(m_branch1);
    aligned_free(m_parent0);
    aligned_free(m_parent1);
}

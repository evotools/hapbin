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

template <bool Binom>
void EHHFinder::calcBranch(const HapMap* hm, HapMap::PrimitiveType* parent, std::size_t parentcount, HapMap::PrimitiveType* branch, std::size_t& branchcount, std::size_t currIndex, double freq, double &probs, std::size_t& singlecount, std::size_t maxBreadth, bool* overflow)
{
    const HapMap::PrimitiveType* mapData = hm->rawData();
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
        
        if (Binom && count <= 1)
        {
            continue;
        }
        if (!Binom && count == 0)
        {
            continue;
        }
        else if (!Binom && count == 1)
        {
            ++singlecount;
            continue;
        }
        else
        {
            if (Binom)
                probs += binom_2(count)*freq;
            else
                probs += (count*freq)*(count*freq);
            for(std::size_t j = 0; j < snpDataSize; ++j)
            {
                branch[bcnt*snpDataSize+j] = parent[i*snpDataSize+j] & mapData[currIndex*snpDataSize+j];
            }
            ++bcnt;
            for(std::size_t j = 0; j < snpDataSize-1; ++j)
            {
                branch[bcnt*snpDataSize+j] = parent[i*snpDataSize+j] & ~mapData[currIndex*snpDataSize+j];
            }
            branch[bcnt*snpDataSize+snpDataSize-1] = (parent[i*snpDataSize+snpDataSize-1] & ~mapData[currIndex*snpDataSize+snpDataSize-1]) & m_maskA;
            ++bcnt;
        }
        if (bcnt > maxBreadth-2)
        {
            *overflow = true;
            return;
        }
    }
    branchcount = bcnt;
}

template <bool Binom>
inline void EHHFinder::calcBranchXPEHH(std::size_t currIndex, std::size_t& singleA, std::size_t& singleB, std::size_t& singleP, bool* overflow)
{
    std::size_t snpDataSize = m_snpDataSizeA + m_snpDataSizeB;
    std::size_t bcnt = 0;
    for(std::size_t i = 0; i < m_parent0count; ++i)
    {
        int countA = 0, countB = 0;
        unsigned long long *leaf = (unsigned long long*) &m_parent0[i*snpDataSize];
        for (std::size_t j = 0; j < m_snpDataSizeULL_A; ++j)
        {
            countA += popcount1(leaf[j]);
        }
        leaf = (unsigned long long*) &m_parent0[i*snpDataSize+ m_snpDataSizeA];
        for (std::size_t j = 0; j < m_snpDataSizeULL_B; ++j)
        {
            countB += popcount1(leaf[j]);
        }
        int count = countA+countB;

        if (Binom && count <= 1)
        {
            continue;
        }
        else if (!Binom && count == 0)
        {
            continue;
        }
        else if (!Binom && count == 1)
        {
            ++singleP;
            singleA+=countA;
            singleB+=countB;
            continue;
        }
        else
        {
            if (Binom)
            {
                m_ehhA += binom_2(countA)*m_freqA;
                m_ehhB += binom_2(countB)*m_freqB;
                m_ehhP += binom_2(count)*m_freqP;
            }
            else
            {
                m_ehhA += (countA*m_freqA)*(countA*m_freqA);
                m_ehhB += (countB*m_freqB)*(countB*m_freqB);
                m_ehhP += (count*m_freqP)*(count*m_freqP);
            }
            //A part of the 1 branch
            for(std::size_t j = 0; j < m_snpDataSizeA; ++j)
            {
                m_branch0[bcnt*snpDataSize+j] = m_parent0[i*snpDataSize+j] & m_hdA[currIndex*m_snpDataSizeA+j];
            }
            //B part of the 1 branch
            for(std::size_t j = 0; j < m_snpDataSizeB; ++j)
            {
                m_branch0[bcnt*snpDataSize+m_snpDataSizeA+j] = m_parent0[i*snpDataSize+m_snpDataSizeA+j] & m_hdB[currIndex*m_snpDataSizeB+j];
            }
            ++bcnt;
            //A part of the 0 branch
            for(std::size_t j = 0; j < m_snpDataSizeA-1; ++j)
            {
                m_branch0[bcnt*snpDataSize+j] = m_parent0[i*snpDataSize+j] & ~m_hdA[currIndex*m_snpDataSizeA+j];
            }
            m_branch0[bcnt*snpDataSize+m_snpDataSizeA-1] = (m_parent0[i*snpDataSize+m_snpDataSizeA-1] & ~m_hdA[(currIndex+1)*m_snpDataSizeA-1]) & m_maskA;
            //B part of the 0 branch
            for(std::size_t j = 0; j < m_snpDataSizeB-1; ++j)
            {
                m_branch0[bcnt*snpDataSize+m_snpDataSizeA+j] = m_parent0[i*snpDataSize+m_snpDataSizeA+j] & ~m_hdB[currIndex*m_snpDataSizeB+j];
            }
            m_branch0[(bcnt+1)*snpDataSize-1] = (m_parent0[(i+1)*snpDataSize-1] & ~m_hdB[(currIndex+1)*m_snpDataSizeB-1]) & m_maskB;
            ++bcnt;
        }
        if (bcnt > m_maxBreadth0-2)
        {
            *overflow = true;
            return;
        }
    }
    m_branch0count = bcnt;
}

template <bool Binom>
void EHHFinder::calcBranchesXPEHH(std::size_t currIndex)
{
    bool overflow;
    bool realloced0 = false;
    std::size_t singleA{}, singleB{}, singleP{};
    do
    {
        overflow = false;
        if (!Binom)
        {
            singleA = 0;
            singleB = 0;
            singleP = 0;
        }
        m_ehhA = 0.0;
        m_ehhB = 0.0;
        m_ehhP = 0.0;
        m_branch0count = 0;
        calcBranchXPEHH<Binom>(currIndex, singleA, singleB, singleP, &overflow);
        if (overflow)
        {
            m_maxBreadth0 += 100;
            aligned_free(m_branch0);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, (m_snpDataSizeA+m_snpDataSizeB)*m_maxBreadth0*sizeof(HapMap::PrimitiveType)));
            realloced0 = true;
            continue;
        }
    } while (false);
    if (realloced0)
    {
        aligned_free(m_parent0);
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, (m_snpDataSizeA+m_snpDataSizeB)*m_maxBreadth0*sizeof(HapMap::PrimitiveType)));
    }
    if (!Binom)
    {
        m_single0count += singleA;
        m_single1count += singleB;
        m_singlePcount += singleP;
    }
    m_parent0count = m_branch0count;
    m_branch0count = 0ULL;
    std::swap(m_parent0, m_branch0);
}

template <bool Binom>
XPEHH EHHFinder::findXPEHH(const HapMap* hmA, const HapMap* hmB, std::size_t focus, std::atomic<unsigned long long>* reachedEnd)
{
    if (focus <= 1 || focus >= hmA->numSnps()-2)
        return XPEHH();
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
    if (!Binom)
    {
        m_single0count = 0ULL;
        m_single1count = 0ULL;
        m_singlePcount = 0ULL;
    }
    unsigned long long locusPysPos = hmA->physicalPosition(focus);
    XPEHH ret;
    ret.index = focus;
    for(std::size_t i = 0; i < m_snpDataSizeA; ++i)
    {
        ret.numA += POPCOUNT(m_hdA[focus*m_snpDataSizeA+i]);
    }
    ret.numNotA = hmA->snpLength() - ret.numA;
    for(std::size_t i = 0; i < m_snpDataSizeB; ++i)
    {
        ret.numB += POPCOUNT(m_hdB[focus*m_snpDataSizeB+i]);
    }
    ret.numNotB = hmA->snpLength() - ret.numB;
    double probASingle, probBSingle, probPSingle;
    double lastEhhA, lastEhhB, lastEhhP;
    if (Binom)
    {
        m_freqA = 1.0/binom_2(hmA->snpLength());
        m_freqB = 1.0/binom_2(hmB->snpLength());
        m_freqP = 1.0/binom_2(hmA->snpLength()+hmB->snpLength());
        lastEhhA = (binom_2(ret.numA)+binom_2(ret.numNotA))*m_freqA;
        lastEhhB = (binom_2(ret.numB)+binom_2(ret.numNotB))*m_freqB;
        lastEhhP = (binom_2(ret.numA+ret.numB)+binom_2(ret.numNotA+ret.numNotB))*m_freqP;
    }
    else
    {
        m_freqA = 1.0/(double)hmA->snpLength();
        m_freqB = 1.0/(double)hmB->snpLength();
        m_freqP = 1.0/(double)(hmA->snpLength()+hmB->snpLength());
        probASingle = m_freqA*m_freqA;
        probBSingle = m_freqB*m_freqB;
        probPSingle = m_freqP*m_freqP;
        double f = ret.numA*m_freqA;
        lastEhhA = f*f+(1.0-f)*(1.0-f);
        f = ret.numB*m_freqB;
        lastEhhB = f*f+(1.0-f)*(1.0-f);
        f = (ret.numA+ret.numB)*m_freqP;
        lastEhhP = f*f+(1.0-f)*(1.0-f);
        m_single0count = 0ULL;
        m_single1count = 0ULL;
    }
    
    if (!Binom)
    {
        m_single0count = 0ULL;
        m_single1count = 0ULL;
        m_singlePcount = 0ULL;
    }
    setInitialXPEHH(focus);
    calcBranchesXPEHH<Binom>(focus-1);
    if (focus > 1)
    {
        for (std::size_t currIndex = focus - 2;; --currIndex)
        { 
            unsigned long long currPhysPos = hmA->physicalPosition(currIndex+1);
            double scale = (double)(m_scale) / (double)(hmA->physicalPosition(currIndex+2) - currPhysPos);
            if (scale > 1)
                scale=1;
            
            calcBranchesXPEHH<Binom>(currIndex);
            
            if (!Binom)
            {
                m_ehhA += probASingle*m_single0count;
                m_ehhB += probBSingle*m_single1count;
                m_ehhP += probPSingle*(m_single0count+m_single1count);
            }

            if (m_ehhP <= m_cutoff - 1e-15)
                break;

            ret.iHH_A1 += (hmA->geneticPosition(currIndex+2)-hmA->geneticPosition(currIndex+1))*(lastEhhA + m_ehhA)*scale*0.5;
            ret.iHH_B1 += (hmA->geneticPosition(currIndex+2)-hmA->geneticPosition(currIndex+1))*(lastEhhB + m_ehhB)*scale*0.5;
            ret.iHH_P1 += (hmA->geneticPosition(currIndex+2)-hmA->geneticPosition(currIndex+1))*(lastEhhP + m_ehhP)*scale*0.5;

            lastEhhA = m_ehhA;
            lastEhhB = m_ehhB;
            lastEhhP = m_ehhP;
            
            if (m_maxExtend != 0 && locusPysPos - currPhysPos > m_maxExtend)
                break;
            if (Binom && m_ehhP == 0)
                break;
            if (!Binom && (m_single0count+m_single1count) == (hmA->snpLength()+hmB->snpLength()))
                break;
            if(currIndex == 0)
            {
                ++(*reachedEnd);
                return XPEHH();
            }
        }
    }
    
    if (Binom)
    {
        lastEhhA = (binom_2(ret.numA)+binom_2(ret.numNotA))*m_freqA;
        lastEhhB = (binom_2(ret.numB)+binom_2(ret.numNotB))*m_freqB;
        lastEhhP = (binom_2(ret.numA+ret.numB)+binom_2(ret.numNotA+ret.numNotB))*m_freqP;
    }
    else
    {
        double f = ret.numA*m_freqA;
        lastEhhA = f*f+(1.0-f)*(1.0-f);
        f = ret.numB*m_freqB;
        lastEhhB = f*f+(1.0-f)*(1.0-f);
        f = (ret.numA+ret.numB)*m_freqP;
        lastEhhP = f*f+(1.0-f)*(1.0-f);
        m_single0count = 0ULL;
        m_single1count = 0ULL;
    }

    setInitialXPEHH(focus);
    calcBranchesXPEHH<Binom>(focus+1);
    for (std::size_t currIndex = focus + 2; currIndex < hmA->numSnps(); ++currIndex)
    {
        unsigned long long currPhysPos = hmA->physicalPosition(currIndex-1);
        double scale = (double)(m_scale) / (double)(currPhysPos - hmA->physicalPosition(currIndex-2));
        if (scale > 1)
            scale=1;
        
        calcBranchesXPEHH<Binom>(currIndex);
        
        if (!Binom)
        {
            m_ehhA += probASingle*m_single0count;
            m_ehhB += probBSingle*m_single1count;
            m_ehhP += probPSingle*(m_single0count+m_single1count);
        }

        if (m_ehhP <= m_cutoff - 1e-15)
            break;

        ret.iHH_A1 += (hmA->geneticPosition(currIndex-1)-hmA->geneticPosition(currIndex-2))*(lastEhhA + m_ehhA)*scale*0.5;
        ret.iHH_B1 += (hmA->geneticPosition(currIndex-1)-hmA->geneticPosition(currIndex-2))*(lastEhhB + m_ehhB)*scale*0.5;
        ret.iHH_P1 += (hmA->geneticPosition(currIndex-1)-hmA->geneticPosition(currIndex-2))*(lastEhhP + m_ehhP)*scale*0.5;

        lastEhhA = m_ehhA;
        lastEhhB = m_ehhB;
        lastEhhP = m_ehhP;
        
        if (m_maxExtend != 0 && currPhysPos - locusPysPos > m_maxExtend)
            break;
        if (!Binom && (m_single0count+m_single1count) == (hmA->snpLength()+hmB->snpLength()))
            break;
        if (Binom && m_ehhP == 0)
            break;
        if (currIndex == hmA->numSnps()-1)
        {
            ++(*reachedEnd);
            return XPEHH();
        }
    }
    return ret;
}

template <bool Binom>
void EHHFinder::calcBranches(const HapMap* hapmap, std::size_t focus, std::size_t currIndex, double freq0,  double freq1, HapStats &stats)
{
    bool overflow;
    bool realloced0 = false, realloced1 = false;
    std::size_t single0{};
    do
    {
        overflow = false;
        if (!Binom)
            single0 = 0;
        stats.probsNot = 0.0;
        calcBranch<Binom>(hapmap, m_parent0, m_parent0count, m_branch0, m_branch0count, currIndex, freq0, stats.probsNot, single0, m_maxBreadth0, &overflow);
        if(overflow)
        {
            m_maxBreadth0 += 100;
            aligned_free(m_branch0);
            m_branch0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth0*sizeof(HapMap::PrimitiveType)));
            realloced0 = true;
            continue;
        }
    } while (false);
    if (realloced0)
    {
        aligned_free(m_parent0);
        m_parent0 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeA*m_maxBreadth0*sizeof(HapMap::PrimitiveType)));
    }
    std::size_t single1{};
    do
    {
        overflow = false;
        if (!Binom)
            single1 = 0;
        stats.probs = 0.0;
        calcBranch<Binom>(hapmap, m_parent1, m_parent1count, m_branch1, m_branch1count, currIndex, freq1, stats.probs, single1, m_maxBreadth1, &overflow);
        if(overflow)
        {
            m_maxBreadth1 += 100;
            aligned_free(m_branch1);
            m_branch1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth1*sizeof(HapMap::PrimitiveType)));
            realloced1 = true;
            continue;
        }
    } while (false);
    if (realloced1)
    {
        aligned_free(m_parent1);
        m_parent1 = reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, m_snpDataSizeB*m_maxBreadth1*sizeof(HapMap::PrimitiveType)));
    }
    m_parent0count = m_branch0count;
    m_parent1count = m_branch1count;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    if (!Binom)
    {
        m_single0count += single0;
        m_single1count += single1;
    }
    std::swap(m_parent0, m_branch0);
    std::swap(m_parent1, m_branch1);
}

template <bool Binom>
EHH EHHFinder::find(const HapMap* hapmap, std::size_t focus, std::atomic<unsigned long long>* reachedEnd, std::atomic<unsigned long long>* outsideMaf, bool ehhsave)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    if (!Binom)
    {
        m_single0count = 0ULL;
        m_single1count = 0ULL;
    }
    m_hdA = hapmap->rawData();
    m_snpDataSizeA = m_snpDataSizeB = hapmap->snpDataSize();

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
    if (!(maxEHH <= 1.0 - m_minMAF && maxEHH >= m_minMAF) && m_minMAF != 0.0)
    {
        ++(*outsideMaf);
        return EHH();
    }
    if (focus < 2 || focus == hapmap->numSnps()-2)
    {
        ++(*reachedEnd);
        return EHH();
    }
    
    double freq0, freq1;
    double probSingle, probNotSingle;
    if (Binom)
    {
        freq0 = 1.0/binom_2(ret.numNot);
        freq1 = 1.0/binom_2(ret.num);
    }
    else
    {
        freq0 = 1.0/(double)ret.numNot;
        freq1 = 1.0/(double)ret.num;
        probSingle = freq1*freq1;
        probNotSingle = freq0*freq0;
    }
    double lastProbs = 1.0, lastProbsNot = 1.0;
    unsigned long long locusPysPos = hapmap->physicalPosition(focus);
    
    setInitial(focus, focus-1);
    
    for (std::size_t currIndex = focus - 2;; --currIndex)
    {
        HapStats stats;
        unsigned long long currPhysPos = hapmap->physicalPosition(currIndex+1);
        double scale = (double)(m_scale) / (double)(hapmap->physicalPosition(currIndex+2) - currPhysPos);
        if (scale > 1)
            scale=1;
        
        calcBranches<Binom>(hapmap, focus, currIndex, freq0, freq1, stats);

        if (!Binom)
        {
            stats.probs +=  probSingle*m_single1count;
            stats.probsNot += probNotSingle*m_single0count;
        }

        if (lastProbs > m_cutoff - 1e-15)
            ret.iHH_1 += (hapmap->geneticPosition(currIndex+2)-hapmap->geneticPosition(currIndex+1))*(lastProbs + stats.probs)*scale*0.5;    
        if (lastProbsNot > m_cutoff - 1e-15)
            ret.iHH_0 += (hapmap->geneticPosition(currIndex+2)-hapmap->geneticPosition(currIndex+1))*(lastProbsNot + stats.probsNot)*scale*0.5;   
        
        lastProbs = stats.probs;
        lastProbsNot = stats.probsNot;
        if (ehhsave)
            ret.upstream.push_back(std::move(stats));
        
        
        if (m_maxExtend != 0 && locusPysPos - currPhysPos > m_maxExtend)
            break;
        if (lastProbs <= m_cutoff - 1e-15 && lastProbsNot <= m_cutoff - 1e-15)
            break;
        if (!Binom && (m_single0count+m_single1count) == hapmap->snpLength())
            break;
        if (currIndex == 0)
        {
            ++(*reachedEnd);
            return EHH();
        }
    }
    
    setInitial(focus,focus+1);
    lastProbs = 1.0, lastProbsNot = 1.0;
    for (std::size_t currIndex = focus + 2; currIndex < hapmap->numSnps(); ++currIndex)
    {
        HapStats stats;
        unsigned long long currPhysPos = hapmap->physicalPosition(currIndex-1);
        double scale = double(m_scale) / double(hapmap->physicalPosition(currIndex-1) - hapmap->physicalPosition(currIndex-2));
        if (scale > 1)
            scale=1;
        
        calcBranches<Binom>(hapmap, focus, currIndex, freq0, freq1, stats);
        
        if (!Binom)
        {
            stats.probs +=  probSingle*m_single1count;
            stats.probsNot += probNotSingle*m_single0count;
        }
        
        if (lastProbs > m_cutoff - 1e-15) {
            ret.iHH_1 += (hapmap->geneticPosition(currIndex-1)-hapmap->geneticPosition(currIndex-2))*(lastProbs + stats.probs)*scale*0.5;
        }
        if (lastProbsNot > m_cutoff - 1e-15) {
            ret.iHH_0 += (hapmap->geneticPosition(currIndex-1)-hapmap->geneticPosition(currIndex-2))*(lastProbsNot + stats.probsNot)*scale*0.5;
        }
        
        lastProbs = stats.probs;
        lastProbsNot = stats.probsNot;
        if (ehhsave)
            ret.downstream.push_back(std::move(stats));
        
        if (lastProbs <= m_cutoff - 1e-15 && lastProbsNot <= m_cutoff - 1e-15 || (m_single0count+m_single1count) == hapmap->snpLength())
            break;
        if (m_maxExtend != 0 && currPhysPos - locusPysPos > m_maxExtend)
            break;
        if (!Binom && currIndex == hapmap->numSnps()-1)
        {
            ++(*reachedEnd);
            return EHH();
        }
    }
    return ret;
}

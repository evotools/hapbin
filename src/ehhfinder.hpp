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

#ifndef LLEHHFINDER_H
#define LLEHHFINDER_H
#include "ehh.hpp"
#include "hapmap.hpp"
#include <atomic>

class EHHFinder
{
public:
    explicit EHHFinder(std::size_t snpDataSizeA, std::size_t snpDataSizeB, std::size_t maxBreadth, double cutoff, double minMAF, double scale, unsigned long long maxExtend);
    EHH find(HapMap* hapmap, std::size_t focus, std::atomic<unsigned long long>* reachedEnd, std::atomic<unsigned long long>* outsideMaf, bool ehhsave = false);
    XPEHH findXPEHH(HapMap* hmA, HapMap *hmB, std::size_t focus, std::atomic<unsigned long long>* reachedEnd);
    ~EHHFinder();
protected:
    inline void calcBranch(HapMap* hm, std::size_t focus, HapMap::PrimitiveType* parent, std::size_t parentcount, HapMap::PrimitiveType* branch, std::size_t& branchcount, std::size_t currLine, double freq, double& probs, std::size_t& singlecount, bool* overflow);
    inline void calcBranchXPEHH(HapMap* hmA, HapMap* hmB, std::size_t currLine, std::size_t& single0, std::size_t& single1, bool* overflow);
    inline void setInitial(std::size_t focus, std::size_t line);
    inline void setInitialXPEHH(std::size_t focus);
    inline void calcBranches(HapMap* hapmap, std::size_t focus, std::size_t currLine, double freq0, double freq1, HapStats& stats);
    inline void calcBranchesXPEHH(std::size_t currLine);
    
    std::size_t m_maxBreadth;
    std::size_t m_bufferSize;
    unsigned long long m_maxExtend;
    HapMap::PrimitiveType *m_parent0;
    HapMap::PrimitiveType *m_parent1;
    HapMap::PrimitiveType *m_branch0;
    HapMap::PrimitiveType *m_branch1;
    HapMap::PrimitiveType m_maskA;
    HapMap::PrimitiveType m_maskB;
    const double m_cutoff;
    const double m_minMAF;
    const double m_scale;
    std::size_t m_maxSnpDataSize;
    std::size_t m_parent0count;
    std::size_t m_parent1count;
    std::size_t m_branch0count;
    std::size_t m_branch1count;
    std::size_t m_single0count;
    std::size_t m_single1count;
    std::size_t m_singlePooledCount;
    double m_freqA;
    double m_freqB;
    double m_freqP;
    double m_ehhA;
    double m_ehhB;
    double m_ehhP;
    std::size_t m_snpDataSizeA;
    std::size_t m_snpDataSizeB;
    std::size_t m_snpDataSizeULL_A;
    std::size_t m_snpDataSizeULL_B;
    HapMap::PrimitiveType *m_hdA;
    HapMap::PrimitiveType *m_hdB;
    HapMap* m_hmA;
    HapMap* m_hmB;
};

#endif // LLEHHFINDER_H

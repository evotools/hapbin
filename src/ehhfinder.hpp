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

#ifndef LLEHHFINDER_H
#define LLEHHFINDER_H
#include "ehh.hpp"
#include "hapmap.hpp"
#include <atomic>

class EHHFinder
{
public:
    explicit EHHFinder(std::size_t snpDataSizeA, std::size_t snpDataSizeB, std::size_t maxBreadth, double cutoff, double minMAF, double scale, unsigned long long maxExtend);
    template <bool Binom>
    EHH find(const HapMap* hapmap, std::size_t focus, std::atomic<unsigned long long>* reachedEnd, std::atomic<unsigned long long>* outsideMaf, bool ehhsave = false);
    template <bool Binom>
    XPEHH findXPEHH(const HapMap* hmA, const HapMap *hmB, std::size_t focus, std::atomic<unsigned long long>* reachedEnd);
    ~EHHFinder();

protected:
    template <bool Binom>
    inline void calcBranch(const HapMap* hm, HapMap::PrimitiveType* parent, std::size_t parentcount, HapMap::PrimitiveType* branch, std::size_t& branchcount, std::size_t currIndex, double freq, double& probs, std::size_t& singlecount, std::size_t maxBreadth, bool* overflow);
    template <bool Binom>
    inline void calcBranchXPEHH(std::size_t currIndex, std::size_t& singleA, std::size_t& singleB, std::size_t& singleP, bool* overflow);
    void setInitial(std::size_t focus, std::size_t index);
    void setInitialXPEHH(std::size_t focus);
    template <bool Binom>
    inline void calcBranches(const HapMap* hapmap, std::size_t focus, std::size_t currIndex, double freq0, double freq1, HapStats& stats);
    template <bool Binom>
    inline void calcBranchesXPEHH(std::size_t currIndex);

    std::size_t m_maxBreadth0;
    std::size_t m_maxBreadth1;
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
    std::size_t m_singlePcount;
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
    const HapMap::PrimitiveType *m_hdA;
    const HapMap::PrimitiveType *m_hdB;
    const HapMap* m_hmA;
    const HapMap* m_hmB;
};

#include "ehhfinder-impl.hpp"

#endif // LLEHHFINDER_H

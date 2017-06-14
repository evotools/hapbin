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

#include "ehhfinder.hpp"
#include "hapmap.hpp"
#include <algorithm>

EHHFinder::EHHFinder(std::size_t snpDataSizeA, std::size_t snpDataSizeB, std::size_t maxBreadth, double cutoff, double minMAF, double scale, unsigned long long maxExtend)
    : m_maxBreadth0(maxBreadth)
    , m_maxBreadth1(maxBreadth)
    , m_bufferSize(snpDataSizeA*maxBreadth)
    , m_maxExtend(maxExtend)
    , m_parent0(reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, (snpDataSizeA+snpDataSizeB)*maxBreadth*sizeof(HapMap::PrimitiveType))))
    , m_parent1(reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, (snpDataSizeA)*maxBreadth*sizeof(HapMap::PrimitiveType))))
    , m_branch0(reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, (snpDataSizeA+snpDataSizeB)*maxBreadth*sizeof(HapMap::PrimitiveType))))
    , m_branch1(reinterpret_cast<HapMap::PrimitiveType*>(aligned_alloc(128, (snpDataSizeA)*maxBreadth*sizeof(HapMap::PrimitiveType))))
    , m_cutoff(cutoff)
    , m_minMAF(minMAF)
    , m_scale(scale)
    , m_maxSnpDataSize(snpDataSizeA)
    , m_freqA{}
    , m_freqB{}
    , m_freqP{}
    , m_ehhA{}
    , m_ehhB{}
    , m_ehhP{}
{}

/**
 * Set initial state. Set m_parent0 to '0' core haplotype positions, m_parent1 to '1' core haplotype positions.
 *
 * calcBranch counts the state of the parent (previous) branch, not the counts of the new (current) level. Therefore, we must advance
 * the calculations by one iteration before starting.
 */
void EHHFinder::setInitial(std::size_t focus, std::size_t index)
{
    m_parent0count = 2ULL;
    m_parent1count = 2ULL;
    m_branch0count = 0ULL;
    m_branch1count = 0ULL;
    m_single0count = 0ULL;
    m_single1count = 0ULL;

    for (std::size_t j = 0; j < m_snpDataSizeA; ++j)
    {
        m_parent0[               j] = ~m_hdA[focus*m_snpDataSizeA+j] &  m_hdA[index*m_snpDataSizeA+j];
        m_parent0[m_snpDataSizeA+j] = (~m_hdA[focus*m_snpDataSizeA+j]) & (~m_hdA[index*m_snpDataSizeA+j]);
        m_parent1[               j] =  m_hdA[focus*m_snpDataSizeA+j] &  m_hdA[index*m_snpDataSizeA+j];
        m_parent1[m_snpDataSizeA+j] =  m_hdA[focus*m_snpDataSizeA+j] & ~m_hdA[index*m_snpDataSizeA+j];
    }
    m_parent0[  m_snpDataSizeA-1] &= m_maskA;
    m_parent0[2*m_snpDataSizeA-1] &= m_maskA;
    m_parent1[2*m_snpDataSizeA-1] &= m_maskA;
}

void EHHFinder::setInitialXPEHH(std::size_t focus)
{
    for (std::size_t i = 0; i < m_snpDataSizeA; ++i)
        m_parent0[i] = ~m_hdA[focus*m_snpDataSizeA+i];
    m_parent0[m_snpDataSizeA-1] &= m_maskA;
    for (std::size_t i = 0; i < m_snpDataSizeB; ++i)
        m_parent0[i+m_snpDataSizeA] = ~m_hdB[focus*m_snpDataSizeB+i];
    m_parent0[m_snpDataSizeA+m_snpDataSizeB-1] &= m_maskB;
    for (std::size_t i = 0; i < m_snpDataSizeA; ++i)
        m_parent0[i+m_snpDataSizeA+m_snpDataSizeB] = m_hdA[focus*m_snpDataSizeA+i];
    for (std::size_t i = 0; i < m_snpDataSizeB; ++i)
        m_parent0[i+2*m_snpDataSizeA+m_snpDataSizeB] = m_hdB[focus*m_snpDataSizeB+i];
}


EHHFinder::~EHHFinder()
{
    aligned_free(m_branch0);
    aligned_free(m_parent0);
    if (m_branch1 != NULL)
        aligned_free(m_branch1);
    if (m_parent1 != NULL)
        aligned_free(m_parent1);
}

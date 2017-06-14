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

#ifndef EHHPAIRFINDER_H
#define EHHPAIRFINDER_H

#include "hapmap.hpp"

struct EHHPair;
class EhhPairFinder
{
public:
    EhhPairFinder(const HapMap* hm, double cutoff, double minMAF, double scale, std::size_t maxBreadth, bool saveEhh = false);
    EHHPair calcEhhPair(std::size_t focus1, std::size_t focus2, bool focus1Conditioned);
    bool passesAf(std::size_t focus1, std::size_t focus2);
    ~EhhPairFinder();

protected:
    inline void setCore(std::size_t focus1, std::size_t focus2);
    inline void createBranches(std::size_t index);
    inline bool processBuffer(std::size_t core, std::size_t index);
    inline void calcBranch(std::size_t core, std::size_t index, bool right);
    inline bool calcBranches(std::size_t index, bool right);
    inline bool branch(std::size_t core, std::size_t index, std::size_t bufferOffset, HapMap::PrimitiveType *parent[4], HapMap::PrimitiveType *branch[4], std::size_t (&branchCount)[4], double freq2[4], double (&ehh)[4], int (&single)[4]);
    inline int countHaps(const HapMap::PrimitiveType* haps);
    inline void calcFreq();
    inline bool checkMaf(std::size_t focus1, std::size_t focus2);
    inline bool reachedCutoff(double ehh);
    inline bool reachedCutoffs();
    inline bool reachedUniqueness(std::size_t core);
    inline bool reachedUniqueness();
    void allocate();
    void reallocBranch();

    const HapMap *m_hm;
    const HapMap::PrimitiveType* m_hd;
    HapMap::PrimitiveType *m_parent[4];
    HapMap::PrimitiveType *m_branch[4];
    HapMap::PrimitiveType *m_parentNot[4];
    HapMap::PrimitiveType *m_branchNot[4];
    double m_ihh[4];
    double m_ihhNot[4];
    double m_sl[4];
    double m_slNot[4];
    double m_ehh[4];
    double m_ehhNot[4];
    double m_lastEhh[4];
    double m_lastEhhNot[4];
    double m_freq2[4];
    double m_freqNot2[4];
    double m_af[4];
    int m_count[4];
    int m_countNot[4];
    int m_single[4];
    int m_singleNot[4];
    bool m_calc[4];
    const double m_cutoff;
    const double m_minMAF;
    const double m_scale;
    std::size_t m_maxBreadth;
    std::size_t m_parentCount[4];
    std::size_t m_branchCount[4];
    std::size_t m_parentNotCount[4];
    std::size_t m_branchNotCount[4];
    std::size_t m_snpDataSize;
    std::size_t m_snpDataSizeULL;
    bool m_saveEhh;
    HapMap::PrimitiveType m_mask;

};

#endif // EHHPAIRFINDER_H

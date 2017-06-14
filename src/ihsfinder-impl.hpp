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
void IHSFinder::runXpehh(const HapMap* mA, const HapMap* mB, std::size_t start, std::size_t end)
{
    #pragma omp parallel shared(mA,mB,start,end)
    {
        EHHFinder finder(mA->snpDataSize(), mB->snpDataSize(), 2000, m_cutoff, m_minMAF, m_scale, m_maxExtend);
        #pragma omp for schedule(dynamic,10)
        for(size_t i = start; i < end; ++i)
        {
            XPEHH xpehh = finder.findXPEHH<Binom>(mA, mB, i, &m_reachedEnd);
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

template <bool Binom>
void IHSFinder::run(const HapMap* map, std::size_t start, std::size_t end)
{
    #pragma omp parallel shared(map, start, end)
    {
        EHHFinder finder(map->snpDataSize(), 0, 2000, m_cutoff, m_minMAF, m_scale, m_maxExtend);
        #pragma omp for schedule(dynamic,10)
        for(size_t i = start; i < end; ++i)
        {
            EHH ehh = finder.find<Binom>(map, i, &m_reachedEnd, &m_outsideMaf);
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

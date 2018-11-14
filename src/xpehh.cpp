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

#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "hapbin.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

void calcXpehhNoMpi(
    const std::string& hapA,
    const std::string& hapB,
    const std::string& map,
    const std::string& outfile,
    double cutoff,
    double minMAF,
    double scale,
    unsigned long long maxExtend,
    int bins,
    bool binom)
{
    HapMap hA, hB;
    if (!hA.loadHap(hapA.c_str()))
    {
        std::cerr << "Error: " << hapA.c_str() << " not found." << std::endl;
        return;
    }
    if (!hB.loadHap(hapB.c_str()))
    {
        std::cerr << "Error: " << hapB.c_str() << " not found." << std::endl;
        return;
    }
#ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    hA.loadMap(map.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    IHSFinder *ihsfinder = new IHSFinder(hA.snpLength() + hB.snpLength(), cutoff, minMAF, scale, maxExtend, bins);
    if (binom)
        ihsfinder->runXpehh<true>(&hA, &hB, 0ULL, hA.numSnps());
    else
        ihsfinder->runXpehh<false>(&hA, &hB, 0ULL, hA.numSnps());

    IHSFinder::LineMap standardized = ihsfinder->normalizeXPEHH();

    auto tend = std::chrono::high_resolution_clock::now();
    auto diff = tend - start;
    std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;

    std::ofstream out(outfile);
    out << "Index\tID\tFreq\tiHH_A1\tiHH_B1\tiHH_P1\tXPEHH\tstd XPEHH" << std::endl;
    for (const auto& it : ihsfinder->unStdXPEHHByLine())
    {
        double freq = (double)(it.second.numA + it.second.numB)/((double) hA.snpLength() + hB.snpLength());
        out << it.first << '\t' << hA.lineToId(it.first) << '\t' << freq << '\t' << it.second.iHH_A1 << '\t' << it.second.iHH_B1 << '\t' << it.second.iHH_P1 << '\t' << it.second.xpehh << '\t' << standardized[it.first] << std::endl;
    }

    std::cout << "# valid loci: " << minMAF << ": " << ihsfinder->unStdXPEHHByLine().size() << std::endl;

    delete ihsfinder;
}

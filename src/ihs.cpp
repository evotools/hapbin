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

void calcIhsNoMpi(
    const std::string& hap,
    const std::string& map,
    const std::string& outfile,
    double cutoff,
    double minMAF,
    double scale,
    unsigned long long maxExtend,
    int bins,
    bool binom)
{
    HapMap hm;
    if (!hm.loadHap(hap.c_str()))
    {
        return;
    }

#ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    hm.loadMap(map.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    IHSFinder *ihsfinder = new IHSFinder(hm.snpLength(), cutoff, minMAF, scale, maxExtend, bins);
    if (binom)
        ihsfinder->run<true>(&hm, 0ULL, hm.numSnps());
    else
        ihsfinder->run<false>(&hm, 0ULL, hm.numSnps());
    IHSFinder::LineMap res = ihsfinder->normalize();

    auto tend = std::chrono::high_resolution_clock::now();
    auto diff = tend - start;
    std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;

    auto unStd = ihsfinder->unStdIHSByLine();
    auto freQs = ihsfinder->freqsByLine();

    std::ofstream out2(outfile);
    out2 << "Index\tID\tFreq\tiHH_0\tiHH_1\tiHS\tStd iHS" << std::endl;

    for (const auto& it : res)
    {
        auto s = unStd[it.first];
        auto q = freQs[it.first];
        out2 << it.first << '\t' << hm.lineToId(it.first) << '\t' << q << '\t' << s.iHH_0 << '\t' << s.iHH_1 << '\t' << s.iHS << "\t" << it.second << std::endl;
    }
    std::cout << "# valid loci: " << res.size() << std::endl;
    std::cout << "# loci with MAF <= " << minMAF << ": " << ihsfinder->numOutsideMaf() << std::endl;
    std::cout << "# loci with NaN result: " << ihsfinder->numNanResults() << std::endl;
    std::cout << "# loci which reached the end of the chromosome: " << ihsfinder->numReachedEnd() << std::endl;
    delete ihsfinder;
}

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

void calcIhsNoMpi(const std::string& hap,
                  const std::string& map,
                  const std::string& outfile,
                  double cutoff,
                  double minMAF,
                  double scale,
                  unsigned long long maxExtend,
                  int bins,
                  bool binom)
{
    HapMap ctcmap;
    if (!ctcmap.loadHap(hap.c_str()))
    {
        return;
    }

#ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    ctcmap.loadMap(map.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    IHSFinder *ihsfinder = new IHSFinder(ctcmap.snpLength(), cutoff, minMAF, scale, maxExtend, bins);
    if (binom)
        ihsfinder->run<true>(&ctcmap, 0ULL, ctcmap.numSnps());
    else
        ihsfinder->run<false>(&ctcmap, 0ULL, ctcmap.numSnps());
    IHSFinder::LineMap res = ihsfinder->normalize();
    
    auto tend = std::chrono::high_resolution_clock::now();
    auto diff = tend - start;
    std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;
    
    auto unStd = ihsfinder->unStdIHSByLine();
    
    /*std::ofstream out(outfile);
    out << "Location\tiHH_0\tiHH_1\tiHS" << std::endl;
    for (const auto& it : unStd)
    {
        out << ctcmap.lineToId(it.first) << '\t' << it.second.iHH_0 << '\t' << it.second.iHH_1 << '\t' << it.second.iHS << std::endl;
    }*/
    std::ofstream out2(outfile);
    out2 << "Location\tiHH_0\tiHH_1\tiHS\tStd iHS" << std::endl;
    
    for (const auto& it : res)
    {
        auto s = unStd[it.first];
        out2 << ctcmap.lineToId(it.first) << '\t' << s.iHH_0 << '\t' << s.iHH_1 << '\t' << s.iHS << "\t" << it.second << std::endl;
    }
    std::cout << "# valid loci: " << res.size() << std::endl;
    std::cout << "# loci with MAF <= " << minMAF << ": " << ihsfinder->numOutsideMaf() << std::endl;
    std::cout << "# loci with NaN result: " << ihsfinder->numNanResults() << std::endl;
    std::cout << "# loci which reached the end of the chromosome: " << ihsfinder->numReachedEnd() << std::endl;
    delete ihsfinder;
}

void calcXpehhNoMpi(const std::string& hapA,
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
    IHSFinder *ihsfinder = new IHSFinder(hA.snpLength(), cutoff, minMAF, scale, maxExtend, bins);
    if (binom)
        ihsfinder->runXpehh<true>(&hA, &hB, 0ULL, hA.numSnps());
    else
        ihsfinder->runXpehh<false>(&hA, &hB, 0ULL, hA.numSnps());
    
    auto tend = std::chrono::high_resolution_clock::now();
    auto diff = tend - start;
    std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;
    
    std::ofstream out(outfile);
    out << "Location\tiHH_A1\tiHH_B1\tiHH_P1\tXPEHH" << std::endl;
    for (const auto& it : ihsfinder->unStdXIHSByLine())
    {
        out << hA.lineToId(it.first) << '\t' << it.second.iHH_A1 << '\t' << it.second.iHH_B1 << '\t' << it.second.iHH_P1 << '\t' << it.second.xpehh << std::endl;
    }

    std::cout << "# valid loci: " << minMAF << ": " << ihsfinder->unStdXIHSByLine().size() << std::endl;

    delete ihsfinder;
}

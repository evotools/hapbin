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

#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "hapbin.hpp"
#include <fstream>

#ifdef _OPENMP
#include <omp.h>
#endif

void calcIhsNoMpi(const HapMap* hm,
                  const std::string& outfile,
                  double cutoff,
                  double minMAF,
                  double scale,
                  unsigned long long maxExtend,
                  int bins,
                  bool binom)
{
#ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    auto start = std::chrono::high_resolution_clock::now();

    IHSFinder *ihsfinder = new IHSFinder(hm->snpLength(), cutoff, minMAF, scale, maxExtend, bins);
    if (binom)
        ihsfinder->run<true>(hm, 0ULL, hm->numSnps());
    else
        ihsfinder->run<false>(hm, 0ULL, hm->numSnps());
    ihsfinder->normalize();
    IHSFinder::IhsInfoMap res = ihsfinder->std();

    auto tend = std::chrono::high_resolution_clock::now();
    auto diff = tend - start;
    std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;

    auto unStd = ihsfinder->unStdByIndex();

    std::ofstream out(outfile);
    out << "Location\tiHH_0\tiHH_1\tiHS\tStd iHS\tSL_0\tSL_1\tnSL\tStd nSL" << std::endl;

    for (const auto& it : res)
    {
        auto s = unStd[it.first];
        out << hm->indexToId(it.first) << '\t';
        out << it.second.iHH_0 << '\t';
        out << it.second.iHH_1 << '\t';
        out << log(it.second.iHH_1/it.second.iHH_0) << "\t";
        out << it.second.iHS << '\t';
        out << it.second.sl_0 << '\t';
        out << it.second.sl_1 << '\t';
        out << log(it.second.sl_1/it.second.sl_1) << '\t';
        out << it.second.nSL << std::endl;
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
                    const std::string& keyA,
                    const std::vector<std::string>& popsA,
                    const std::string& keyB,
                    const std::vector<std::string>& popsB,
                    const std::string& outfile,
                    double cutoff,
                    double minMAF,
                    double scale,
                    unsigned long long maxExtend,
                    int bins,
                    bool binom)
{
    HapMap mA, mB;
    PopKey *pkA = NULL, *pkB = NULL;
    if (keyA.size() > 0 && popsA.size() > 0)
        pkA = new PopKey(keyA, popsA);
    if (keyB.size() > 0 && popsB.size() > 0)
        pkB = new PopKey(keyB, popsB);
    if (!mA.loadHap(hapA.c_str(), false, pkA))
    {
        std::cerr << "Error: " << hapA.c_str() << " not found." << std::endl;
        return;
    }
    if (!mB.loadHap(hapB.c_str(), false, pkB))
    {
        std::cerr << "Error: " << hapB.c_str() << " not found." << std::endl;
        return;
    }
#ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    mA.loadMap(map.c_str());
    auto start = std::chrono::high_resolution_clock::now();
    IHSFinder *ihsfinder = new IHSFinder(mA.snpLength(), cutoff, minMAF, scale, maxExtend, bins);
    if (binom)
        ihsfinder->runXpehh<true>(&mA, &mB, 0ULL, mA.numSnps());
    else
        ihsfinder->runXpehh<false>(&mA, &mB, 0ULL, mA.numSnps());

    auto tend = std::chrono::high_resolution_clock::now();
    auto diff = tend - start;
    std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;


    std::ofstream out(outfile);
    out << "Location\tiHH_A1\tiHH_B1\tiHH_P1\tXPEHH\tsl_A1\tsl_B1\tsl_P1\tXPnSl" << std::endl;
    for (const auto& it : ihsfinder->unStdXIHSByIndex())
    {
        out << mA.indexToId(it.first) << '\t' << it.second.iHH_A1 << '\t' << it.second.iHH_B1 << '\t' << it.second.iHH_P1 << '\t' << it.second.xpehh << '\t' << it.second.sl_A1 << '\t' << it.second.sl_B1 << '\t' << it.second.sl_P1 << '\t' << log(it.second.sl_A1/it.second.sl_B1) << std::endl;
    }

    std::cout << "# valid loci: " << minMAF << ": " << ihsfinder->unStdByIndex().size() << std::endl;

    delete ihsfinder;
}

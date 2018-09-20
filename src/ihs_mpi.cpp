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

#include "config.h"
#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "hapbin.hpp"
#include "ehh.hpp"

#if MPI_FOUND
#include "mpirpc/manager.hpp"
#include "mpirpc/parameterstream.hpp"

ParameterStream& operator<<(ParameterStream& out, const IhsScore& info)
{
    out << info.iHS << info.iHH_0 << info.iHH_1 << info.freq;
    return out;
}

ParameterStream& operator>>(ParameterStream& in, IhsScore& info)
{
    in >> info.iHS >> info.iHH_0 >> info.iHH_1 >> info.freq;
    return in;
}
#endif

#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif


void calcIhsMpi(
    const std::string& hapfile,
    const std::string& mapfile,
    const std::string& outfile,
    double cutoff,
    double minMAF,
    double scale,
    unsigned long long maxExtend,
    int binFactor,
    bool binom)
{
#if MPI_FOUND
    std::cout << "Calculating iHS using MPI." << std::endl;
    HapMap hap;
    if (!hap.loadHap(hapfile.c_str()))
    {
        return;
    }
    std::cout << "Loaded " << hap.numSnps() << " snps." << std::endl;
    std::cout << "Haplotype count: " << hap.snpLength() << " " << maxExtend << std::endl;
    hap.loadMap(mapfile.c_str());
    IHSFinder *ihsfinder = new IHSFinder(hap.snpLength(), cutoff, minMAF, scale, maxExtend, binFactor);
    mpirpc::Manager *manager = new mpirpc::Manager();
    int procsToGo = manager->numProcs();
    std::cout << "Processes: " << procsToGo << std::endl;
#ifdef _OPENMP
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    mpirpc::FunctionHandle done = manager->registerLambda([&]() {
        --procsToGo;
    });
    manager->registerType<IHSFinder>();
    manager->registerFunction(&IHSFinder::addData);
    if (manager->rank() == 0)
    {
        manager->registerObject(ihsfinder);
    }
    /**
     * MPI_Issend() in Manager::registerObject does not necessarily notify other processes that a send is ready before the barrier.
     * Therefore, we must loop until the sends and recieves are complete.
     */
    while (manager->getObjectsOfType<IHSFinder>().size() < 1 || manager->queueSize() > 0)
    {
        manager->checkMessages();
    }
    manager->barrier();
    mpirpc::ObjectWrapperBase* mainihsfinder = *(manager->getObjectsOfType<IHSFinder>().cbegin());
    mpirpc::FunctionHandle runEHH =  manager->registerLambda([&](std::size_t start, std::size_t end) {
        std::cout << "Computing EHH for " << (end-start) << " lines on rank " << manager->rank() << std::endl;
        if (binom)
            ihsfinder->run<true>(&hap, start, end);
        else
            ihsfinder->run<false>(&hap, start, end);
        IHSFinder::LineMap fbl = ihsfinder->freqsByLine();
        manager->invokeFunction(mainihsfinder, &IHSFinder::addData, false, fbl, ihsfinder->unStdIHSByLine(), ihsfinder->unStdIHSByFreq(), ihsfinder->numReachedEnd(), ihsfinder->numOutsideMaf(), ihsfinder->numNanResults());
        manager->invokeFunction(0, done);
    });
    manager->barrier();

    auto start = std::chrono::high_resolution_clock::now();

    if (manager->rank() == 0)
    {
        std::size_t numSnps = hap.numSnps();
        std::size_t snpsPerRank = numSnps/manager->numProcs();
        std::size_t pos = 0ULL;
        for (int i = 1; i < manager->numProcs(); ++i)
        {
            std::size_t start = snpsPerRank*i;
            std::size_t end = snpsPerRank*(i+1);
            if (i == manager->numProcs()-1)
            {
                end = numSnps;
            }
            manager->invokeFunction(i, runEHH, start, end);
            manager->checkMessages();
        }
        while (manager->queueSize() > 0) {
            manager->checkMessages();
        }
        std::size_t end = snpsPerRank;
        if (manager->numProcs() == 1)
            end = numSnps;
        if (binom)
            ihsfinder->run<true>(&hap, 0UL, end);
        else
            ihsfinder->run<false>(&hap, 0UL, end);
        --procsToGo;
    }

    while(manager->checkMessages())
    {
        if (procsToGo == 0)
            manager->shutdown();
    }

    if (manager->rank() == 0)
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto diff = end - start;
        std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;

        IHSFinder::LineMap res = ihsfinder->normalize();

        auto unStd = ihsfinder->unStdIHSByLine();

        /*std::ofstream out(outfile);
        out << "Location\tiHH_0\tiHH_1\tiHS" << std::endl;
        for (const auto& it : unStd)
        {
            out << hap.lineToId(it.first) << '\t' << it.second.iHH_0 << '\t' << it.second.iHH_1 << '\t' << it.second.iHS << std::endl;
        }*/
        std::ofstream out2(outfile);
        out2 << "Location\tFreq\tiHH_0\tiHH_1\tiHS\tStd iHS" << std::endl;

        for (const auto& it : res)
        {
            auto s = unStd[it.first];
            out2 << hap.lineToId(it.first) << '\t' << s.freq << '\t' << s.iHH_0 << '\t' << s.iHH_1 << '\t' << s.iHS << "\t" << it.second << std::endl;
        }
        std::cout << "# valid loci: " << res.size() << std::endl;
        std::cout << "# loci with MAF <= " << minMAF << ": " << ihsfinder->numOutsideMaf() << std::endl;
        std::cout << "# loci with NaN result: " << ihsfinder->numNanResults() << std::endl;
        std::cout << "# loci which reached the end of the chromosome: " << ihsfinder->numReachedEnd() << std::endl;
    }

    delete ihsfinder;
    delete manager;
#else
    std::cout << "MPI support not enabled!" << std::endl;
#endif
}

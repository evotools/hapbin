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
#include "ehh.hpp"

#if MPI_FOUND
#include "mpirpc/manager.hpp"
#include "mpirpc/parameterstream.hpp"

ParameterStream& operator<<(ParameterStream& out, const XPEHH& info)
{
    out << info.index << info.xpehh << info.numA << info.numB << info.numNotA << info.numNotB << info.iHH_A1 << info.iHH_B1 << info.iHH_P1;
    return out;
}

ParameterStream& operator>>(ParameterStream& in, XPEHH& info)
{
    in >> info.index >> info.xpehh >> info.numA >> info.numB >> info.numNotA >> info.numNotB >>  info.iHH_A1 >> info.iHH_B1 >> info.iHH_P1;
    return in;
}
#endif

#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

void calcXpehhMpi(
    const std::string& hapA,
    const std::string& hapB,
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
    std::cout << "Calculating XPEHH using MPI." << std::endl;
    HapMap mA, mB;
    if (!mA.loadHap(hapA.c_str()))
    {
        return;
    }
    if (!mB.loadHap(hapB.c_str()))
    {
        return;
    }
    std::cout << "Loaded " << mA.numSnps() << " snps for population A." << std::endl;
    std::cout << "Loaded " << mB.numSnps() << " snps for population B." << std::endl;
    std::cout << "Population A haplotype count: " << mA.snpLength() << std::endl;
    std::cout << "Population B haplotype count: " << mB.snpLength() << std::endl;
    mA.loadMap(mapfile.c_str());
    IHSFinder *ihsfinder = new IHSFinder(mA.snpLength() + mB.snpLength(), cutoff, minMAF, scale, maxExtend, binFactor);
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
    manager->registerFunction(&IHSFinder::addXData);
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
    mpirpc::FunctionHandle runXPEHH =  manager->registerLambda([&](std::size_t start, std::size_t end) {
        std::cout << "Computing XPEHH for " << (end-start) << " lines on rank " << manager->rank() << std::endl;
        if (binom)
            ihsfinder->runXpehh<true>(&mA, &mB, start, end);
        else
            ihsfinder->runXpehh<false>(&mA, &mB, start, end);
        IHSFinder::LineMap fbl = ihsfinder->freqsByLine();
        manager->invokeFunction(mainihsfinder, &IHSFinder::addXData, false, fbl, ihsfinder->unStdXPEHHByLine(), ihsfinder->unStdIHSByFreq(), ihsfinder->numReachedEnd(), ihsfinder->numOutsideMaf(), ihsfinder->numNanResults());
        manager->invokeFunction(0, done);
    });
    manager->barrier();

    auto start = std::chrono::high_resolution_clock::now();

    if (manager->rank() == 0)
    {

        std::size_t numSnps = mA.numSnps();
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
            manager->invokeFunction(i, runXPEHH, start, end);
            manager->checkMessages();
        }
        while (manager->queueSize() > 0) {
            manager->checkMessages();
        }
        std::size_t end = snpsPerRank;
        if (manager->numProcs() == 1)
            end = numSnps;
        if (binom)
            ihsfinder->runXpehh<true>(&mA, &mB, 0ULL, end);
        else
            ihsfinder->runXpehh<false>(&mA, &mB, 0ULL, end);
        --procsToGo;
    }

    while(manager->checkMessages())
    {
        if (procsToGo == 0)
            manager->shutdown();
    }

    if (manager->rank() == 0)
    {
        IHSFinder::LineMap standardized = ihsfinder->normalizeXPEHH();

        auto tend = std::chrono::high_resolution_clock::now();
        auto diff = tend - start;
        std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;

        std::ofstream out(outfile);
        out << "Index\tID\tFreq\tiHH_A1\tiHH_B1\tiHH_P1\tXPEHH\tstd XPEHH" << std::endl;
        for (const auto& it : ihsfinder->unStdXPEHHByLine())
        {
            double freq = (double)(it.second.numA + it.second.numB)/((double) mA.snpLength() + mB.snpLength());
            out << it.first << '\t' << mA.lineToId(it.first) << '\t' << freq << '\t' << it.second.iHH_A1 << '\t' << it.second.iHH_B1 << '\t' << it.second.iHH_P1 << '\t' << it.second.xpehh << '\t' << standardized[it.first] << std::endl;
        }
        std::cout << "# valid loci: " << ihsfinder->unStdXPEHHByLine().size() << std::endl;
    }
    delete ihsfinder;
    delete manager;
#else
    std::cout << "MPI support not enabled!" << std::endl;
#endif
}

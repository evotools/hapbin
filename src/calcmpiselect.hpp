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

#include <mpirpc/manager.hpp>
#include <mpirpc/parameterstream.hpp>

#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

void calcIhsMpi(const std::string& hapfile, const std::string& mapfile, const std::string& outfile, double cutoff, double minMAF, double scale, double binFactor)
{
    std::cout << "Calculating iHS using MPI." << std::endl;
    HapMap hap;
    if (!hap.loadHap(hapfile.c_str()))
    {
        return;
    }
    std::cout << "Loaded " << hap.numSnps() << " snps." << std::endl;
    std::cout << "Haplotype count: " << hap.snpLength() << std::endl;
    hap.loadMap(mapfile.c_str());
    IHSFinder *ihsfinder = new IHSFinder(hap.snpLength(), cutoff, minMAF, scale, binFactor);
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
        ihsfinder->run(&hap, start, end);
        IHSFinder::LineMap fbl = ihsfinder->freqsByLine();
        manager->invokeFunction(mainihsfinder, &IHSFinder::addData, false, fbl, ihsfinder->unStdIHSByLine(), ihsfinder->unStdIHSByFreq());
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
        ihsfinder->run(&hap, 0UL, end);
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
    
        std::ofstream out(outfile);
        for (const auto& it : ihsfinder->unStdIHSByLine())
        {
            out << hap.lineToId(it.first) << " " << it.second << std::endl;
        }
        std::ofstream out2(outfile+".std");
        for (const auto& it : res)
        {
            out2 << hap.lineToId(it.first) << " " << it.second << std::endl;
        }
        std::cout << "# loci with MAF >= " << minMAF << ": " << res.size() << std::endl;
    }
    
    delete ihsfinder;
    delete manager;
}

void calcXpehhMpi(const std::string& hapA, const std::string& hapB, const std::string& mapfile, const std::string& outfile, double cutoff, double minMAF, double scale, double binFactor)
{
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
    IHSFinder *ihsfinder = new IHSFinder(mA.snpLength(), cutoff, minMAF, scale, binFactor);
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
    mpirpc::FunctionHandle runXPEHH =  manager->registerLambda([&](std::size_t start, std::size_t end) {
        std::cout << "Computing XPEHH for " << (end-start) << " lines on rank " << manager->rank() << std::endl;
        ihsfinder->runXpehh(&mA, &mB, start, end);
        IHSFinder::LineMap fbl = ihsfinder->freqsByLine();
        manager->invokeFunction(mainihsfinder, &IHSFinder::addData, false, fbl, ihsfinder->unStdIHSByLine(), ihsfinder->unStdIHSByFreq());
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
        ihsfinder->runXpehh(&mA, &mB, 0ULL, end);
        --procsToGo;
    }
    
    while(manager->checkMessages())
    {
        if (procsToGo == 0)
            manager->shutdown();
    }
    
    if (manager->rank() == 0)
    {
        auto tend = std::chrono::high_resolution_clock::now();
        auto diff = tend - start;
        std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;
        
        std::ofstream out(outfile);
        for (const auto& it : ihsfinder->unStdIHSByLine())
        {
            out << mA.lineToId(it.first) << " " << it.second << std::endl;
        }
    }
    delete ihsfinder;
    delete manager;
}

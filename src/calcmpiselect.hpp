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
#include "ihspairfinder.hpp"
#include "hapbin.hpp"
#include "ehh.hpp"
#include "ihspairjob.hpp"

#include "mpirpc/manager.hpp"
#include "mpirpc/parameterstream.hpp"
#include "mpirpc/mpitype.hpp"

#include <ctime>
#include <fstream>
#include <list>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace impl
{
    template<template<typename,typename,typename...> class C, typename T, typename U, typename...O>
    inline void bufferSize(const C<T,U,O...>& v, std::size_t& size)
    {
        size += v.size();
    }

    template<template<typename,typename,typename...> class C, typename T, typename U, typename...O>
    inline void fillBuffer(const C<T,U,O...>& v, U* buf, std::size_t& index)
    {
        for (const auto& p : v)
            buf[index++] = p.second;
    }

    template<template<typename,typename,typename...> class C, typename T, typename U, typename... O>
    inline void fillMap(C<T,U,O...>& v, const U* buf, std::size_t&index)
    {
        for (auto& p : v)
            p.second = buf[index++];
    }
}

template<template<typename,typename,typename...> class... C, typename T, typename U, typename... O>
void reduceMaps(MPI_Op op, MPI_Comm comm, C<T,U,O...>&... in)
{
    std::size_t sumBufferSize = 0ULL, index = 0ULL;
    mpirpc::Passer{(impl::bufferSize(in, sumBufferSize), 0)...};
    U *sumBuffer    = new U[sumBufferSize];
    U *sumBufferRes = new U[sumBufferSize];
    mpirpc::Passer{(impl::fillBuffer(in, sumBuffer, index), 0)...};
    MPI_Allreduce(sumBuffer, sumBufferRes, sumBufferSize, mpiType<U>(), op, comm);
    index = 0ULL;
    mpirpc::Passer{(impl::fillMap(in,sumBufferRes,index), 0)...};
    delete[] sumBuffer;
    delete[] sumBufferRes;
}

ParameterStream& operator<<(ParameterStream& out, const IhsScore& info)
{
    out << info.iHS << info.iHH_0 << info.iHH_1 << info.nSL << info.sl_0 << info.sl_1;
    return out;
}

ParameterStream& operator>>(ParameterStream& in, IhsScore& info)
{
    in >> info.iHS >> info.iHH_0 >> info.iHH_1 >> info.nSL >> info.sl_0 >> info.sl_1;
    return in;
}

ParameterStream operator<<(ParameterStream& out, const StdScorePair& s)
{
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        out << s.freq[core];
    }
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        out << s.score[core];
    }
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        out << s.score_0[core];
    }
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        out << s.score_1[core];
    }
    return out;
}

ParameterStream operator>>(ParameterStream& in, StdScorePair& s)
{
    double freq[4], score[4], score_0[4], score_1[4];
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        in >> freq[core];
    }
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        in >> score[core];
    }
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        in >> score_0[core];
    }
    for(std::size_t core = 0b00; core <= 0b11; ++core)
    {
        in >> score_1[core];
    }
    return in;
}

ParameterStream& operator<<(ParameterStream& out, const XPEHH& info)
{
    out << info.index << info.xpehh << info.numA << info.numB << info.numNotA << info.numNotB << info.iHH_A1 << info.iHH_B1 << info.iHH_P1 << info.sl_A1 << info.sl_B1 << info.sl_P1;
    return out;
}

ParameterStream& operator>>(ParameterStream& in, XPEHH& info)
{
    in >> info.index >> info.xpehh >> info.numA >> info.numB >> info.numNotA >> info.numNotB >>  info.iHH_A1 >> info.iHH_B1 >> info.iHH_P1 >> info.sl_A1 >> info.sl_B1 >> info.sl_P1;
    return in;
}

void calcIhsMpi(const HapMap* hap,
                const std::string& outfile,
                double cutoff,
                double minMAF,
                double scale,
                unsigned long long maxExtend,
		        int binFactor,
		        bool binom)
{
    std::cout << "Calculating iHS using MPI." << std::endl;
    std::cout << "Loaded " << hap->numSnps() << " snps." << std::endl;
    std::cout << "Haplotype count: " << hap->snpLength() << " " << maxExtend << std::endl;
    IHSFinder *ihsfinder = new IHSFinder(hap->snpLength(), cutoff, minMAF, scale, maxExtend, binFactor);
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
    manager->sync();
    mpirpc::ObjectWrapperBase* mainihsfinder = *(manager->getObjectsOfType<IHSFinder>().cbegin());
    mpirpc::FunctionHandle runEHH =  manager->registerLambda([&](std::size_t start, std::size_t end) {
        std::cout << "Computing EHH for " << (end-start) << " lines on rank " << manager->rank() << std::endl;
        if (binom)
            ihsfinder->run<true>(hap, start, end);
        else
            ihsfinder->run<false>(hap, start, end);
        IHSFinder::IndexMap fbl = ihsfinder->freqsByIndex();
        manager->invokeFunction(mainihsfinder, &IHSFinder::addData, false, fbl, ihsfinder->unStdByIndex(), ihsfinder->unStdIHSByFreq(), ihsfinder->unStdNSLByFreq(), ihsfinder->numReachedEnd(), ihsfinder->numOutsideMaf(), ihsfinder->numNanResults());
        manager->invokeFunction(0, done);
    });
    manager->sync();

    auto start = std::chrono::high_resolution_clock::now();

    if (manager->rank() == 0)
    {
        std::size_t numSnps = hap->numSnps();
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
            ihsfinder->run<true>(hap, 0UL, end);
        else
            ihsfinder->run<false>(hap, 0UL, end);
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

        ihsfinder->normalize();
        IHSFinder::IhsInfoMap std = ihsfinder->std();

        std::ofstream out(outfile);
        out << "Location\tiHH_0\tiHH_1\tiHS\tStd iHS\tSL_0\tSL_1\tnSL\tStd nSL" << std::endl;

        for (const auto& it : std)
        {
            double unStdIHS = log(it.second.iHH_1/it.second.iHH_0);
            double unStdNSL = log(it.second.sl_1/it.second.sl_0);
            out << hap->indexToId(it.first) << '\t';
            out << it.second.iHH_0 << '\t';
            out << it.second.iHH_1 << '\t';
            out << unStdIHS << "\t";
            out << it.second.iHS << '\t';
            out << it.second.sl_0 << '\t';
            out << it.second.sl_1 << '\t';
            out << unStdNSL << '\t';
            out << it.second.nSL << std::endl;
        }
        std::cout << "# valid loci: " << std.size() << std::endl;
        std::cout << "# loci with MAF <= " << minMAF << ": " << ihsfinder->numOutsideMaf() << std::endl;
        std::cout << "# loci with NaN result: " << ihsfinder->numNanResults() << std::endl;
        std::cout << "# loci which reached the end of the chromosome: " << ihsfinder->numReachedEnd() << std::endl;
    }

    delete ihsfinder;
    delete manager;
}

void calcXpehhMpi(const std::string& hapA,
                  const std::string& hapB,
                  const std::string& mapfile,
                  const std::string& keyA,
                  const std::vector<std::string>& popsA,
                  const std::string& keyB,
                  const std::vector<std::string>& popsB,
                  const std::string& outfile,
                  double cutoff,
                  double minMAF,
                  double scale,
                  unsigned long long maxExtend,
		          int binFactor,
		          bool binom)
{
    std::cout << "Calculating XPEHH using MPI." << std::endl;
    HapMap mA, mB;
    PopKey *pkA = NULL, *pkB = NULL;
    if (keyA.size() > 0 && popsA.size() > 0)
        pkA = new PopKey(keyA, popsA);
    if (keyB.size() > 0 && popsB.size() > 0)
        pkB = new PopKey(keyB, popsB);
    if (!mA.loadHap(hapA.c_str(), false, pkA))
    {
        return;
    }
    if (!mB.loadHap(hapB.c_str(), false, pkB))
    {
        return;
    }
    std::cout << "Loaded " << mA.numSnps() << " snps for population A." << std::endl;
    std::cout << "Loaded " << mB.numSnps() << " snps for population B." << std::endl;
    std::cout << "Population A haplotype count: " << mA.snpLength() << std::endl;
    std::cout << "Population B haplotype count: " << mB.snpLength() << std::endl;
    mA.loadMap(mapfile.c_str());
    IHSFinder *ihsfinder = new IHSFinder(mA.snpLength(), cutoff, minMAF, scale, maxExtend, binFactor);
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
    manager->sync();
    mpirpc::ObjectWrapperBase* mainihsfinder = *(manager->getObjectsOfType<IHSFinder>().cbegin());
    mpirpc::FunctionHandle runXPEHH =  manager->registerLambda([&](std::size_t start, std::size_t end) {
        std::cout << "Computing XPEHH for " << (end-start) << " lines on rank " << manager->rank() << std::endl;
        if (binom)
            ihsfinder->runXpehh<true>(&mA, &mB, start, end);
        else
            ihsfinder->runXpehh<false>(&mA, &mB, start, end);
        IHSFinder::IndexMap fbl = ihsfinder->freqsByIndex();
        manager->invokeFunction(mainihsfinder,
                                &IHSFinder::addXData,
                                false,
                                fbl,
                                ihsfinder->unStdXIHSByIndex(),
                                ihsfinder->unStdIHSByFreq(),
                                ihsfinder->numReachedEnd(),
                                ihsfinder->numOutsideMaf(),
                                ihsfinder->numNanResults());
        manager->invokeFunction(0, done);
    });
    manager->sync();

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
        auto tend = std::chrono::high_resolution_clock::now();
        auto diff = tend - start;
        std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;

        std::ofstream out(outfile);
        out << "Location\tiHH_A1\tiHH_B1\tiHH_P1\tXPEHH\tsl_A1\tsl_B1\tsl_P1\tXPnSl" << std::endl;
        for (const auto& it : ihsfinder->unStdXIHSByIndex())
        {
            out << mA.indexToId(it.first) << '\t' << it.second.iHH_A1 << '\t' << it.second.iHH_B1 << '\t' << it.second.iHH_P1 << '\t' << it.second.xpehh << '\t' << it.second.sl_A1 << '\t' << it.second.sl_B1 << '\t' << it.second.sl_P1 << '\t' << log(it.second.sl_A1/it.second.sl_B1) << std::endl;
        }
        std::cout << "# valid loci: " << ihsfinder->unStdByIndex().size() << std::endl;
    }
    delete ihsfinder;
    delete manager;
}

void calcIhs2Mpi(const std::string& hapFile,
                 const std::string& mapFile,
                 bool filter, const PopKey* pk,
                 const std::string& key,
                 const std::vector<std::string>& pops,
                 double cutoff,
                 double minMaf,
                 double scale,
                 unsigned long long window,
                 std::size_t maxBreadth,
                 double sig,
                 std::size_t jobSize,
                 unsigned int bins,
                 unsigned int windowbins,
                 const std::string& outfile,
                 unsigned long long xsize,
                 unsigned long long ysize,
                 unsigned long long startPos = 0ULL,
                 unsigned long long endPos = 0ULL)
{
    mpirpc::Manager manager;
    HapMap *hm = nullptr;
    IhsPairFinder *finder = nullptr;
    bool loaded = false;
    std::list<std::pair<std::size_t, std::size_t>> *jobPool;
    std::mutex jobPoolMutex;
    bool finished = false;
    int slavesRunning = manager.numProcs()-1;
    std::size_t totalJobs = 0;
    std::string outprefix = outfile;
    struct timespec waitTime, t2;
    waitTime.tv_sec  = 0;
    waitTime.tv_nsec = 10000000L;
    std::cout << "Processes: " << manager.numProcs() << std::endl;
#ifdef _OPENMP
    omp_set_nested(true);
    std::cout << "Threads: " << omp_get_max_threads() << std::endl;
#endif
    std::size_t jobStart;
    jobPool = new std::list<std::pair<std::size_t, std::size_t>>();
    std::map<std::pair<std::size_t,std::size_t>,StdScorePair> significant_ihs_pairs;
    std::map<std::pair<std::size_t,std::size_t>,StdScorePair> significant_nsl_pairs;
    std::vector<IhsPairJob*> ownedJobs;
    if (manager.rank() == 0)
    {
        hm = new HapMap();
        hm->loadHap(hapFile.c_str(), filter, pk);
        hm->loadMap(mapFile.c_str());
        std::cout << "Loaded " << hm->numSnps() << " SNPs." << std::endl;
        finder = new IhsPairFinder(hm, cutoff, minMaf, scale, window, maxBreadth);
        //std::cout << "Total number of pairs: " << finder->numPairs() << std::endl;
        loaded = true;
        if (endPos == 0)
            endPos = hm->numSnps();
        for(jobStart = startPos; jobStart < endPos-jobSize; jobStart+=jobSize)
        {
            jobPool->push_back({jobStart, jobStart+jobSize});
        }
        jobPool->push_back({jobStart, endPos});
        totalJobs = jobPool->size();
    }
    mpirpc::FunctionHandle createFinder = manager.registerLambda([&](const std::string& hapFile, const std::string& mapFile, bool filter, const std::string& key, const std::vector<std::string>& pops, double cutoff, double minMaf, double scale, unsigned long long window, std::size_t maxBreadth){
        hm = new HapMap();
        PopKey *popkey = NULL;
        if (key.size() > 0 && pops.size() > 0)
            popkey = new PopKey(key, pops);
        hm->loadHap(hapFile.c_str(), filter, popkey);
        hm->loadMap(mapFile.c_str());
        finder = new IhsPairFinder(hm, cutoff, minMaf, scale, window, maxBreadth);
        loaded = true;
    });
    //Set variables which may not have been set (MPI implementations may only pass program args to rank 0)
    mpirpc::FunctionHandle setVars = manager.registerLambda([&](int b, const std::string& prefix){
        bins = b;
        outprefix = prefix;
    });
    mpirpc::FunctionHandle done = manager.registerLambda([&]() {
        --slavesRunning;
    });
    mpirpc::FunctionHandle getJob = manager.registerLambda([&](int remoteRank) -> std::pair<std::size_t,std::size_t> {
        jobPoolMutex.lock();
        if(jobPool->size() > 0)
        {
            auto ret = jobPool->back();
            std::cout << "Processing job " << ret.first << "-" << ret.second << " on rank " << remoteRank << ": " << jobPool->size() << "/" << totalJobs << " jobs left." << std::endl;
            jobPool->pop_back();
            jobPoolMutex.unlock();
            return ret;
        }
        else
        {
            jobPoolMutex.unlock();
            return {0,0};
        }
    });
    mpirpc::FunctionHandle send_significant_ihs = manager.registerLambda([&](const std::map<std::pair<std::size_t,std::size_t>,StdScorePair> m){
        significant_ihs_pairs.insert(m.cbegin(),m.cend());
    });
    mpirpc::FunctionHandle send_significant_nsl = manager.registerLambda([&](const std::map<std::pair<std::size_t,std::size_t>,StdScorePair> m){
        significant_nsl_pairs.insert(m.cbegin(),m.cend());
    });
    manager.sync();
    unsigned long long total_pairs = 0;
    auto start = std::chrono::high_resolution_clock::now();
    if (manager.rank() == 0)
    {
        for(int i = 1; i < manager.numProcs(); ++i)
        {
            manager.invokeFunction(i, createFinder, hapFile, mapFile, filter, key, pops, cutoff, minMaf, scale, window, maxBreadth);
            manager.invokeFunction(i, setVars, bins, outprefix);
        }
        #pragma omp parallel num_threads(2)
        {
            #pragma omp master
            {
                while(manager.checkMessages())
                {
                    if (slavesRunning == 0)
                        break;
                    nanosleep(&waitTime, &t2);
                }
            } //no implied barrier
            #pragma omp single //Necessary for 1 MPI procs
            {
                while (jobPool->size() > 0)
                {
                    std::pair<std::size_t, std::size_t> job{};
                    jobPoolMutex.lock();
                    if (jobPool->size() > 0)
                    {
                        job = jobPool->back();
                        std::cout << "Processing job " << job.first << "-" << job.second << " on rank 0: " << jobPool->size() << "/" << totalJobs << " jobs left." << std::endl;
                        jobPool->pop_back();
                    }
                    jobPoolMutex.unlock();
                    if (job.first != 0 || job.second != 0)
                    {
                        IhsPairJob *pre = new IhsPairJob();
                        if (pre->load(outprefix + "." + std::to_string(job.first) + "-" + std::to_string(job.second) + ".dat") == 0)
                        {
                            total_pairs += pre->size();
                            std::cout << "Reusing results for " << job.first << "-" << job.second << ", which has " << pre->size() << " pairs. # previously calculated pairs, total: " << total_pairs << std::endl;
                            //finder->add(pre.results());
                            ownedJobs.push_back(pre);
                            continue;
                        }
                        else
                        {
                            delete pre;
                        }
                        auto start_time = std::chrono::high_resolution_clock::now();
                        IhsPairJob *res = finder->calcRange(job.first, job.second);
                        auto end_time = std::chrono::high_resolution_clock::now();
                        std::cout << "Job " << job.first << "-" << job.second << " took " << std::chrono::duration<double, std::milli>(end_time-start_time).count() << "ms." << std::endl;
                        res->save(outprefix + "." + std::to_string(job.first) + "-" + std::to_string(job.second) + ".dat");
                        ownedJobs.push_back(res);
                        //delete res;
                    }
                }
            }
        } //implied barrier
    }
    else
    {
        while(manager.checkMessages())
        {
            if (loaded)
            {
                auto job = manager.invokeFunction<std::pair<std::size_t,std::size_t>>(0, getJob, true, manager.rank());
                if (job.first == 0 && job.second == 0)
                {
                    finished = true;
                    break;
                }
                else
                {
                    IhsPairJob *pre = new IhsPairJob();
                    if (pre->load(outprefix + "." + std::to_string(job.first) + "-" + std::to_string(job.second) + ".dat") == 0)
                    {
                        total_pairs += pre->size();
                        std::cout << "Reusing results for " << job.first << "-" << job.second << ", which has " << pre->size() << " pairs. # previously calculated pairs, total: " << total_pairs << std::endl;
                        //finder->add(pre.results());
                        ownedJobs.push_back(pre);
                        continue;
                    }
                    else
                    {
                        delete pre;
                    }
                    auto start_time = std::chrono::high_resolution_clock::now();
                    IhsPairJob *res = finder->calcRange(job.first, job.second);
                    res->save(outprefix + "." + std::to_string(job.first) + "-" + std::to_string(job.second) + ".dat");
                    auto end_time = std::chrono::high_resolution_clock::now();
                    std::cout << "Job " << job.first << "-" << job.second << " took " << std::chrono::duration<double, std::milli>(end_time-start_time).count() << "ms." << std::endl;
                    ownedJobs.push_back(pre);
                    //delete res;
                }
            }
        }
        std::cout << "MPI process " << manager.rank() << " finished calculations." << std::endl;
        manager.invokeFunction(0, done);
    }
    manager.sync();
    if (manager.rank() == 0)
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto diff = end - start;
        std::cout << "Calculations took " << std::chrono::duration<double, std::milli>(diff).count() << "ms" << std::endl;
    }
    unsigned long long total_pairs_all = 0;
    MPI_Reduce(&total_pairs,&total_pairs_all,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,manager.comm());
    if (manager.rank() == 0)
        std::cout << "Total number of pairs across all ranks: " << total_pairs_all << std::endl;

    //assert(finder->results().size() == 0);
    /*size_t numJobs = hm->numSnps()/jobSize;
    size_t numProcs = manager.numProcs();
    size_t jobsPerProc = numJobs/numProcs;
    size_t startSnp = jobsPerProc*jobSize*manager.rank();
    size_t endSnp = jobsPerProc*jobSize*(manager.rank()+1);*/
    /*std::cout << "Jobs per proc: " << jobsPerProc << " Start SNP: " << startSnp << " End SNP: " << endSnp << "rank: " << manager.rank() << std::endl;
    if (manager.rank() != manager.numProcs() - 1)
        for(jobStart = startSnp; jobStart < endSnp; jobStart+=jobSize)
        {
            IhsPairJob f;
            f.load(outprefix + "." + std::to_string(jobStart) + "-" + std::to_string(jobStart+jobSize) + ".dat");
            finder->add(f.results());
            std::cout << "Rank " << manager.rank() << " loaded " << outprefix << "." << std::to_string(jobStart) << "-" << std::to_string(jobStart+jobSize) << ".dat" << std::endl;
        }
    else
    {
        endSnp = hm->numSnps();
        for(jobStart = startSnp; jobStart < hm->numSnps()-jobSize; jobStart+=jobSize)
        {
            IhsPairJob f;
            f.load(outprefix + "." + std::to_string(jobStart) + "-" + std::to_string(jobStart+jobSize) + ".dat");
            finder->add(f.results());
            std::cout << "Last rank, " << manager.rank() << ", loaded " << outprefix << "." << std::to_string(jobStart) << "-" << std::to_string(jobStart+jobSize) << ".dat" << std::endl;
        }
        IhsPairJob f;
        f.load(outprefix + "." + std::to_string(jobStart) + "-" + std::to_string(hm->numSnps()) + ".dat");
        finder->add(f.results());
        std::cout << "Last rank, " << manager.rank() << ", loaded tail " << outprefix << "." << std::to_string(jobStart) << "-" << std::to_string(jobStart+jobSize) << ".dat" << std::endl;
    }*/

    //auto unStd = finder->results();

    //std::cout << "Finished loading results" << std::endl;

    /*std::ofstream out(outprefix + "." + std::to_string(manager.rank()) + ".unStd");
    for (const IhsPair& p : unStd)
    {
        out << p.focus[0] << "\t" << p.focus[1] << "\t";
        for (std::size_t core = 0b00; core <= 0b11; ++core)
        {
            out << p.af[core] << "\t" << p.ihh[core] << "\t" << p.ihhNot[core] << "\t" << log(p.ihh[core]/p.ihhNot[core]) << "\t";
        }
        out << std::endl;
    }*/
    double *grid = new double[xsize*ysize]{};
    std::map<std::pair<int,double>,double> unStdIhsSumByFreq;
    std::map<std::pair<int,double>,double> unStdIhsSumSqByFreq;
    std::map<std::pair<int,double>,unsigned long long> unStdIhsCountByFreq;
    std::map<std::pair<int,double>,double> unStdIhsMeanByFreq;
    std::map<std::pair<int,double>,double> unStdIhsStdDevByFreq;
    std::map<std::pair<int,double>,double> unStdNSlSumByFreq;
    std::map<std::pair<int,double>,double> unStdNSlSumSqByFreq;
    std::map<std::pair<int,double>,unsigned long long> unStdNSlCountByFreq;
    std::map<std::pair<int,double>,double> unStdNSlMeanByFreq;
    std::map<std::pair<int,double>,double> unStdNSlStdDevByFreq;
    //Ensure each rank has all the bins
    for (unsigned int i = 0; i < bins; ++i)
    {
        for (unsigned int j = 0; j < windowbins; ++j)
        {
            unStdIhsSumByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdIhsSumSqByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdIhsMeanByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdIhsStdDevByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdIhsCountByFreq[std::pair<int,double>(j,i/(double)bins)] = 0;

            unStdNSlSumByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdNSlSumSqByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdNSlMeanByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdNSlStdDevByFreq[std::pair<int,double>(j,i/(double)bins)] = 0.0;
            unStdNSlCountByFreq[std::pair<int,double>(j,i/(double)bins)] = 0;
        }
    }
    for (const auto& job : ownedJobs)
    {
        for (const auto& p : job->results())
        {
            int windowbin = ((hm->physicalPosition(p.focus[1])-hm->physicalPosition(p.focus[0]))/(double) window)*windowbins;
            for(std::size_t core = 0b00; core <= 0b11; ++core)
            {
                double unstdIhs = log(p.ihh[core]/p.ihhNot[core]);
                if (std::isinf(unstdIhs) || std::isnan(unstdIhs))
                    continue;
                assert(p.af[core] <= 1.0 && p.af[core] >= 0);
                double freq = ((int) (bins*p.af[core]))/(double)bins;
                unStdIhsSumByFreq[{windowbin,freq}] += unstdIhs;
                ++unStdIhsCountByFreq[{windowbin,freq}];
            }
            for(std::size_t core = 0b00; core <= 0b11; ++core)
            {
                double unstdNSl = log(p.sl[core]/p.slNot[core]);
                if (std::isinf(unstdNSl) || std::isnan(unstdNSl))
                    continue;
                assert(p.af[core] <= 1.0 && p.af[core] >= 0);
                double freq = ((int) (bins*p.af[core]))/(double)bins;
                unStdNSlSumByFreq[{windowbin,freq}] += unstdNSl;
                ++unStdNSlCountByFreq[{windowbin,freq}];
            }
        }
    }

    std::cout << "Finished sums on rank " << manager.rank() << std::endl;
    manager.sync();

    reduceMaps(MPI_SUM, manager.comm(), unStdIhsSumByFreq, unStdNSlSumByFreq);
    reduceMaps(MPI_SUM, manager.comm(), unStdIhsCountByFreq, unStdNSlCountByFreq);

    manager.sync();

    if (manager.rank() == 0)
        std::cout << "Finished all sums" << std::endl;

    for (auto& p : unStdIhsMeanByFreq)
        p.second = unStdIhsSumByFreq[p.first]/(double)unStdIhsCountByFreq[p.first];
    for (auto& p : unStdNSlMeanByFreq)
        p.second = unStdNSlSumByFreq[p.first]/(double)unStdNSlCountByFreq[p.first];

    for (const auto& job : ownedJobs)
    {
        for (const auto& p : job->results())
        {
            int windowbin = ((hm->physicalPosition(p.focus[1])-hm->physicalPosition(p.focus[0]))/(double) window)*windowbins;
            for(std::size_t core = 0b00; core <= 0b11; ++core)
            {
                double unstdIhs = log(p.ihh[core]/p.ihhNot[core]);
                if (std::isinf(unstdIhs) || std::isnan(unstdIhs))
                    continue;
                assert(p.af[core] <= 1.0 && p.af[core] >= 0);
                double freq = ((int) (bins*p.af[core]))/(double)bins;
                unStdIhsSumSqByFreq[{windowbin,freq}] += std::pow(unstdIhs-unStdIhsMeanByFreq[{windowbin,freq}], 2);
            }
            for(std::size_t core = 0b00; core <= 0b11; ++core)
            {
                double unstdNSl = log(p.sl[core]/p.slNot[core]);
                if (std::isinf(unstdNSl) || std::isnan(unstdNSl))
                   continue;
                assert(p.af[core] <= 1.0 && p.af[core] >= 0);
                double freq = ((int) (bins*p.af[core]))/(double)bins;
                unStdNSlSumSqByFreq[{windowbin,freq}] += std::pow(unstdNSl-unStdNSlMeanByFreq[{windowbin,freq}], 2);
            }
        }
    }

    std::cout << "Finished sum sqs on rank " << manager.rank() << std::endl;
    manager.sync();

    reduceMaps(MPI_SUM, manager.comm(), unStdIhsSumSqByFreq, unStdNSlSumSqByFreq);

    manager.sync();
    if (manager.rank() == 0)
        std::cout << "Finished sum sqs on all ranks" << std::endl;

    for (const auto& p : unStdIhsSumSqByFreq)
        unStdIhsStdDevByFreq[p.first] = std::sqrt(p.second/(double)unStdIhsCountByFreq[p.first]);
    for (const auto& p : unStdNSlSumSqByFreq)
        unStdNSlStdDevByFreq[p.first] = std::sqrt(p.second/(double)unStdNSlCountByFreq[p.first]);

    std::cout << "Writing stats on rank" << manager.rank() << std::endl;

    if (manager.rank() == 0)
    {
        std::ofstream statsout(outprefix + ".stats");
        statsout << "Window Bin\tFrequency Bin\tCount\tMean\tStd. Dev." << std::endl;
        auto it1 = unStdIhsCountByFreq.cbegin();
        auto it2 = unStdIhsMeanByFreq.cbegin();
        auto it3 = unStdIhsStdDevByFreq.cbegin();
        for(; it1 != unStdIhsCountByFreq.cend() && it2 != unStdIhsMeanByFreq.cend() && it3 != unStdIhsStdDevByFreq.cend(); ++it1, ++it2, ++it3)
        {
            unsigned long long count      = it1->second;
            unsigned long long window_bin = it1->first.first;
            double freq_bin               = it1->first.second;
            double mean     = it2->second;
            double std_dev  = it3->second;
            statsout << (window/windowbins)*window_bin << "\t" << freq_bin << "\t" << count << "\t" << mean << "\t" << std_dev << std::endl;
        }
        statsout.close();
    }
    for (const auto& job : ownedJobs)
    {
        std::ofstream out2(outprefix + "." + std::to_string(job->start()) + "-" + std::to_string(job->end()) + ".std.dat", std::ios::out | std::ios::binary);
        uint64_t count = job->size();
        out2.write(reinterpret_cast<const char*>(&IHS2_STD_MNUM), sizeof(uint64_t));
        out2.write(reinterpret_cast<const char*>(&count), sizeof(uint64_t));
        for (const auto& p : job->results())
        {
            uint64_t f0 = p.focus[0], f1 = p.focus[1];
            int windowbin = ((hm->physicalPosition(p.focus[1])-hm->physicalPosition(p.focus[0]))/(double) window)*windowbins;
            std::size_t gridpos_x = ((p.focus[0])/(double) (hm->numSnps()))*xsize;
            std::size_t gridpos_y = ((p.focus[1])/(double) (hm->numSnps()))*ysize;
            //std::size_t gridpos_y = ((hm->physicalPosition(p.focus[1])-hm->physicalPosition(p.focus[0]))/(double) window)*ysize;
            out2.write(reinterpret_cast<const char*>(&f0), sizeof(uint64_t));
            out2.write(reinterpret_cast<const char*>(&f1), sizeof(uint64_t));
            double tmp_high_ihs[4];
            double tmp_high_nsl[4];
            bool has_high_ihs = false;
            bool has_high_nsl = false;
            for (std::size_t core = 0b00; core <= 0b11; ++core)
            {
                double freq = ((int) (bins*p.af[core]))/(double)bins;
                double unstdIhs = log(p.ihh[core]/p.ihhNot[core]);
                double stdIhs = (unstdIhs - unStdIhsMeanByFreq.at({windowbin,freq}))/unStdIhsStdDevByFreq.at({windowbin,freq});
                double unstdNSl = log(p.sl[core]/p.slNot[core]);
                double stdNSl = (unstdNSl - unStdNSlMeanByFreq.at({windowbin,freq}))/unStdNSlStdDevByFreq.at({windowbin,freq});
                if (!std::isnan(stdIhs) && !std::isinf(stdIhs))
                    grid[gridpos_y*ysize+gridpos_x] += stdIhs*stdIhs*stdIhs;
                tmp_high_ihs[core] = stdIhs;
                tmp_high_nsl[core] = stdNSl;
                if (stdIhs > sig || stdIhs < -sig)
                    has_high_ihs = true;
                if (stdNSl > sig || stdNSl < -sig)
                    has_high_nsl = true;
                out2.write(reinterpret_cast<const char*>(&(p.af[core])), sizeof(double));
                out2.write(reinterpret_cast<const char*>(&stdIhs), sizeof(double));
                out2.write(reinterpret_cast<const char*>(&stdNSl), sizeof(double));
            }
            if (has_high_ihs)
            {
                significant_ihs_pairs.emplace(std::make_pair(std::pair<std::size_t,std::size_t>(p.focus[0],p.focus[1]),StdScorePair(p.af, tmp_high_ihs, p.ihhNot, p.ihh)));
            }
            if (has_high_nsl)
            {
                significant_nsl_pairs.emplace(std::make_pair(std::pair<std::size_t,std::size_t>(p.focus[0],p.focus[1]),StdScorePair(p.af, tmp_high_nsl, p.slNot, p.sl)));
            }
        }
        out2.close();
    }

    std::cout << "Finished writing stats on rank " << manager.rank() << std::endl;

    manager.sync();

    /*MPI_Datatype local_vec;
    MPI_Type_vector(ysize,local_xsize,local_xsize,MPI_DOUBLE,&local_vec);
    MPI_Type_commit(&local_vec);

    MPI_Datatype vec;
    MPI_Type_vector(ysize,local_xsize,xsize,MPI_DOUBLE,&vec);
    MPI_Type_commit(&vec);

    int *displs  = new int[manager.numProcs()]{0};
    int *rcounts = new int[manager.numProcs()]{1};
    for (std::size_t i = 0; i < manager.numProcs(); ++i)
    {
        displs[i] = i*local_xsize;
    }
    double *grid = nullptr;
    if (manager.rank() == 0)
        grid = new double[ysize*local_xsize]{};
    MPI_Gatherv(local_grid,1,local_vec,grid,rcounts,displs,vec,0,manager.comm());*/

    double* grid2 = nullptr;
    if (manager.rank() == 0)
        grid2 = new double[xsize*ysize]{};
    MPI_Reduce(grid,grid2,xsize*ysize,MPI_DOUBLE,MPI_SUM,0,manager.comm());

    manager.sync();

    delete[] grid;

    if (manager.rank() == 0)
    {
        double max_score = 0.0;
        double min_score = 0.0;
        for (std::size_t i = 0; i < ysize; ++i)
        {
            for(std::size_t j = 0; j < xsize; ++j)
            {
                if (grid2[i*ysize+j] > max_score)
                    max_score = grid2[i*ysize+j];
                if (grid2[i*ysize+j] < min_score)
                    min_score = grid2[i*ysize+j];
            }
        }
        std::cout << "max_score: " << max_score << std::endl;
        std::cout << "min_score: " << min_score << std::endl;

        FILE *ppm_file = fopen((outprefix + ".ppm").c_str(),"wb");
        (void) fprintf(ppm_file, "P6\n%d %d\n65535\n", (int) xsize, (int) ysize);
        for (std::size_t i = 0; i < ysize; ++i)
        {
            for(std::size_t j = 0; j < xsize; ++j)
            {
                unsigned short color[3]{};
                if (grid2[i*ysize+j] > 0)
                {
                    color[0] = 65535;
                    color[1] = 65535*(1.0-grid2[i*ysize+j]/max_score);
                    color[2] = 65535*(1.0-grid2[i*ysize+j]/max_score);
                    (void) fwrite(color,sizeof(short),3,ppm_file);
                }
                else
                {
                    color[0] = 65535*(1.0-grid2[i*ysize+j]/min_score);
                    color[1] = 65535*(1.0-grid2[i*ysize+j]/min_score);
                    color[2] = 65535;
                    (void) fwrite(color,sizeof(short),3,ppm_file);
                }
                //std::cout << grid2[i*ysize+j] << " " << color[0] << " " << color[1] << " " << color[2] << std::endl;
            }
        }
        (void) fclose(ppm_file);
        std::cout << "Finished writing image" << std::endl;
    }

    //Process nSL in parallel if we can
    int nsl_output_rank = (manager.numProcs() > 1) ? 1 : 0;

    if (manager.rank() > 0)
    {
        manager.invokeFunction(0,send_significant_ihs,significant_ihs_pairs);
    }
    if (manager.rank() != nsl_output_rank)
    {
        manager.invokeFunction(nsl_output_rank,send_significant_nsl,significant_nsl_pairs);
    }
    manager.sync();
    if (manager.rank() == 0)
    {
        std::ofstream out_s_ihs(outprefix + ".sig.ihs");
        out_s_ihs << "focus0\tfocus1\t"
                     "freq00\tiHH_0,00\tiHH_1,00\tStd. iHS_00\t"
                     "freq01\tiHH_0,01\tiHH_1,01\tStd. iHS_01\t"
                     "freq10\tiHH_0,10\tiHH_1,10\tStd. iHS_10\t"
                     "freq11\tiHH_0,11\tiHH_1,11\tStd. iHS_11" << std::endl;
        for (const auto& kv : significant_ihs_pairs)
        {
            out_s_ihs << hm->indexToId(kv.first.first) << "\t" << hm->indexToId(kv.first.second) << "\t";
            for (std::size_t core = 0b00; core <= 0b10; ++core)
            {
                out_s_ihs << kv.second.freq[core] << "\t";
                out_s_ihs << kv.second.score_0[core] << "\t";
                out_s_ihs << kv.second.score_1[core] << "\t";
                out_s_ihs << kv.second.score[core] << "\t";
            }
            out_s_ihs << kv.second.freq[0b11] << "\t";
            out_s_ihs << kv.second.score_0[0b11] << "\t";
            out_s_ihs << kv.second.score_1[0b11] << "\t";
            out_s_ihs << kv.second.score[0b11] << std::endl;
        }
        out_s_ihs.close();
    }
    if (manager.rank() == nsl_output_rank)
    {
        std::ofstream out_s_nsl(outprefix + ".sig.nsl");
        out_s_nsl << "focus0\tfocus1\t"
                     "freq00\tSL_0,00\tSL_1,00\tStd. nSL_00\t"
                     "freq01\tSL_0,01\tSL_1,01\tStd. nSL_01\t"
                     "freq10\tSL_0,10\tSL_1,10\tStd. nSL_10\t"
                     "freq11\tSL_0,11\tSL_1,11\tStd. nSL_11" << std::endl;
        for (const auto& kv : significant_nsl_pairs)
        {
            out_s_nsl << hm->indexToId(kv.first.first) << "\t" << hm->indexToId(kv.first.second) << "\t";
            for (std::size_t core = 0b00; core <= 0b10; ++core)
            {
                out_s_nsl << kv.second.freq[core] << "\t";
                out_s_nsl << kv.second.score_0[core] << "\t";
                out_s_nsl << kv.second.score_1[core] << "\t";
                out_s_nsl << kv.second.score[core] << "\t";
            }
            out_s_nsl << kv.second.freq[0b11] << "\t";
            out_s_nsl << kv.second.score_0[0b11] << "\t";
            out_s_nsl << kv.second.score_1[0b11] << "\t";
            out_s_nsl << kv.second.score[0b11] << std::endl;
        }
        out_s_nsl.close();
    }
    manager.sync();
    if (manager.rank() == 0)
        std::cout << "Cleaning up" << std::endl;
    delete finder;
    if (manager.rank() == 0)
        delete[] grid2;
    for (auto job : ownedJobs)
        delete job;
    delete hm;
    if (manager.rank() == 0)
        delete jobPool;
    if (manager.rank() == 0)
        std::cout << "Done" << std::endl;
}

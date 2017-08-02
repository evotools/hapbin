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

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "config.h"
#include "argparse.hpp"
#include "ehhpairfinder.hpp"
#include "ihspairfinder.hpp"
#include "ihspairjob.hpp"

#include <functional>
#include <cstdlib>
#include <vector>
#include <utility>
#include <mutex>
#include <fstream>

int main(int argc, char** argv)
{
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<std::string> hap('d', "hap", "Hap file", false, true, "");
    Argument<std::string> map('m', "map", "Map file", false, true, "");
    Argument<std::string> list('l', "list", "List of pairs to calculate. Format: <ID1> <ID2>", false, true, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('b', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<unsigned long> scale('s', "scale", "Gap scale parameter in pb, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<unsigned long long> window('w', "window", "Maximum gap in bp between two core SNP pairs.", false, false, 500000);
    Argument<unsigned long> bufferSize('u', "buffer", "Maximum starting buffer size in SNPs. Will grow if necessary.", false, false, 2000);
    Argument<unsigned int> numBins('i', "bins", "Number of frequency bins to use for iHS standardization", false, false, 50);
    Argument<std::string> outfile('o', "out", "Output file", false, false, "out");
    Argument<bool> filternonpoly('f', "filter", "Filter out non-polymorphic loci", false, false);
    Argument<std::string> pops('p', "pop", "Filter by population", true, false,"");
    Argument<std::string> key('k', "key", "Population key", false, false, "");
    ArgParse argparse({&help, &hap, &map, &list, &outfile, &cutoff, &minMAF, &scale, &window, &bufferSize, &numBins, &filternonpoly, &pops, &key},
                      "Usage: ihs2binsub --map input.map --hap input.hap --list list.txt [--ascii] [--out outfile]");
    if (!argparse.parseArguments(argc, argv))
        return 1;

    std::size_t numSnps = HapMap::querySnpLength(hap.value().c_str());
    std::cout << "Chromosomes per SNP: " << numSnps << std::endl;
    HapMap hm;
    PopKey *popkey = NULL;
    if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
    {
        std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
        return 3;
    }
    if (pops.wasFound() && key.wasFound())
        popkey = new PopKey(key.value(), pops.values());
    hm.loadHap(hap.value().c_str(), filternonpoly.value(), popkey);
    hm.loadMap(map.value().c_str());
    std::ifstream list_file(list.value());
    std::string line;
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    while (std::getline(list_file, line)) {
        std::vector<std::string> split = splitString(line, ' ');
        if (split.size() == 1)
            split = splitString(line, '\t');
        if (split.size() == 3) {
            pairs.emplace_back(std::make_pair(hm.idToIndex(split[0].c_str()), hm.idToIndex(split[1].c_str())));
        }
    }
    std::cout << "# Pairs: " << pairs.size() << std::endl;
    IhsPairJob *job = new IhsPairJob();
    std::mutex mutex;
    #pragma omp parallel shared(pairs)
    {
        EhhPairFinder pf(&hm, cutoff.value(), minMAF.value(), scale.value(), bufferSize.value());
        #pragma omp for schedule(dynamic,1)
        for(auto i = pairs.begin(); i < pairs.end(); ++i)
        {
            EHHPair res = pf.calcEhhPair(i->first, i->second, false);
            if (res.focus[0] != 0UL && res.focus[1] != 0UL)
            {
                mutex.lock();
                job->add(res);
                mutex.unlock();
            }
        }
    }
    job->saveAscii(outfile.value(),hm);

    delete popkey;
    return 0;
}



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
#include "calcselect.hpp"
#include "config.h"
#include "argparse.hpp"
#include "ehhpairfinder.hpp"
#include "ihspairfinder.hpp"

#if MPI_FOUND
#include <mpi.h>
#endif
#include <functional>
#include <cstdlib>

int main(int argc, char** argv)
{
#if MPI_FOUND
    MPI_Init(&argc, &argv);
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<const char*> hap('d', "hap", "Hap file", false, true, "");
    Argument<const char*> map('m', "map", "Map file", false, true, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('b', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<double> sig('x', "sig", "Create filtered iHS/nSL output for pairs with standardized iHS/nSL at least this value (default: 4.5)", false, false, 4.0);
    Argument<unsigned long> scale('s', "scale", "Gap scale parameter in pb, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<unsigned long long> window('w', "window", "Maximum gap in bp between two core SNP pairs.", false, false, 1000000);
    Argument<unsigned long> bufferSize('u', "buffer", "Maximum starting buffer size in SNPs. Will grow if necessary.", false, false, 2000);
    Argument<unsigned long> jobSize('j', "job-size", "Number of first focus SNPs to take at a time for each MPI process.", false, false, 100);
    Argument<unsigned int> numBins('i', "bins", "Number of frequency bins to use for iHS standardization", false, false, 50);
    Argument<unsigned int> windowBins('x', "window-bins", "Number of window bins to use for iHS standardization", false, false, 50);
    Argument<std::string> outfile('o', "out", "Output file", false, false, "out");
    Argument<bool> filternonpoly('f', "filter", "Filter out non-polymorphic loci", false, false);
    Argument<std::string> pops('p', "pop", "Filter by population", true, false,"");
    Argument<std::string> key('k', "key", "Population key", false, false, "");
    Argument<unsigned long long> startPos('y', "start-index", "Index to start pair jobs", false, false, 0);
    Argument<unsigned long long> endPos('z', "end-index", "Index to end pair jobs (0 for end)", false, false, 0);
    Argument<unsigned long long> xsize(ArgumentBase::NO_SHORT_OPT, "x-image-size", "Image x dimension size", false, false, 4096);
    Argument<unsigned long long> ysize(ArgumentBase::NO_SHORT_OPT, "y-image-size", "Image y dimension size", false, false, 4096);
    ArgParse argparse({&help, &hap, &map, &outfile, &cutoff, &minMAF, &sig, &scale, &window, &bufferSize, &jobSize, &numBins, &windowBins, &filternonpoly, &pops, &key, &startPos, &endPos, &xsize, &ysize},
                      "Usage: ihs2bin --map input.map --hap input.hap [--ascii] [--out outfile]");
    if (!argparse.parseArguments(argc, argv))
        return 1;

    std::size_t numSnps = HapMap::querySnpLength(hap.value());
    std::cout << "Chromosomes per SNP: " << numSnps << std::endl;

    PopKey *popkey = NULL;
    if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
    {
        std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
        return 3;
    }
    if (pops.wasFound() && key.wasFound())
        popkey = new PopKey(key.value(), pops.values());

    calcIhs2Mpi(hap.value(), map.value(), filternonpoly.value(), popkey, key.value(), pops.values(), cutoff.value(), minMAF.value(), scale.value(), window.value(),
                bufferSize.value(), sig.value(), jobSize.value(), numBins.value(), windowBins.value(), outfile.value(), xsize.value(), ysize.value(), startPos.value(), endPos.value());

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
#else
    std::cout << "ERROR: ihs2bin must be built with MPI and MPIRPC!" << std::endl;
#endif
    return 0;
}



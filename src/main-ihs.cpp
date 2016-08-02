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

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ihsfinder.hpp"
#include "calcselect.hpp"
#include "config.h"
#include "argparse.hpp"

#if MPI_FOUND
#include <mpi.h>
#endif
#include <functional>
#include <cstdlib>

int main(int argc, char** argv)
{
    int ret = 0;
    std::size_t numSnps;
#if MPI_FOUND
    MPI_Init(&argc, &argv);
#endif
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<bool> version('v', "version", "Version information", true, false);
    Argument<std::string> hap('d', "hap", "Hap file", false, false, "");
    Argument<std::string> map('m', "map", "Map file", false, false, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('f', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<int> binfac('b', "bin", "Number of frequency bins for iHS normalization (default: 50)", false, false, 50);
    Argument<unsigned long long> scale('s', "scale", "Gap scale parameter in bp, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<bool> binom('a', "binom", "Use binomial coefficients rather than frequency squared for EHH", true, false);
    Argument<std::string> outfile('o', "out", "Output file", false, false, "out.txt");
    ArgParse argparse({&help, &version, &hap, &map, &outfile, &cutoff, &minMAF, &scale, &binfac, &binom}, "Usage: ihsbin --map input.map --hap input.hap [--ascii] [--out outfile]");
    if (!argparse.parseArguments(argc, argv)) 
    {
        ret = 1;
        goto out;
    }
    if (help.value())
    {
        argparse.showHelp();
        goto out;
    }
    else if (version.value())
    {
        argparse.showVersion();
        goto out;
    }
    else if (!hap.wasFound() || !map.wasFound())
    {
        std::cout << "Please specify --hap and --map." << std::endl;
        ret = 2;
        goto out;
    }
    
    numSnps = HapMap::querySnpLength(hap.value().c_str());
    std::cout << "Chromosomes per SNP: " << numSnps << std::endl;
    
    calcIhs(hap.value(), map.value(), outfile.value(), cutoff.value(), minMAF.value(), (double) scale.value(), binfac.value(), binom.value());

out:
#if MPI_FOUND
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
#endif
    return ret;
}



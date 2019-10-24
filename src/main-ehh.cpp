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

#include <iostream>

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ehh.hpp"
#include "argparse.hpp"
#include "ehhfinder.hpp"
#include <atomic>

int main(int argc, char** argv)
{
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<bool> version('v', "version", "Version information", true, false);
    Argument<const char*> hap('d', "hap", "Hap file", false, false, "");
    Argument<const char*> map('m', "map", "Map file", false, false, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('b', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<unsigned long long> scale('s', "scale", "Gap scale parameter in bp, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<bool> binom('a', "binom", "Use binomial coefficients rather than frequency squared for EHH", true, false);
    Argument<unsigned long long> maxExtend('e', "max-extend", "Maximum distance in bp to traverse when calculating EHH (default: 0 (disabled))", false, false, 0);
    Argument<const char*> locus('l', "locus", "Locus", false, false, 0);
    ArgParse argparse({&help, &version, &hap, &map, &locus, &cutoff, &minMAF, &scale, &maxExtend, &binom}, "Usage: ehhbin --map input.map --hap input.hap --locus id");
    if (!argparse.parseArguments(argc, argv)) 
    {
        return 3;
    }
    if (help.value())
    {
        argparse.showHelp();
        return 0;
    }
    else if (version.value())
    {
        argparse.showVersion();
        return 0;
    }
    else if (!hap.wasFound() || !map.wasFound() || !locus.wasFound())
    {
        std::cout << "Please specify --hap, --map, and --locus." << std::endl;
        return 4;
    }
    using HapMapType = HapMap;
    HapMapType hmap;
    if (!hmap.loadHap(hap.value()))
    {
        return 1;
    }
    hmap.loadMap(map.value());
    EHH e;
    std::size_t l = hmap.idToLine(locus.value());
    if (l == std::numeric_limits<std::size_t>::max())
    {
        std::cerr << "no locus with the id: " << locus.value() << std::endl;
        return 2;
    }
    std::atomic<unsigned long long> reachedEnd{};
    std::atomic<unsigned long long> outsideMaf{};
    EHHFinder finder(hmap.snpDataSize(), 0, 1000, cutoff.value(), minMAF.value(), (double) scale.value(), maxExtend.value());
    if (binom.value())
        e = finder.find<true>(&hmap, l, &reachedEnd, &outsideMaf, true);
    else
        e = finder.find<false>(&hmap, l, &reachedEnd, &outsideMaf, true);
    e.printEHH(&hmap);
    std::cout << "IHH_0: " << e.iHH_0 << " IHH_1: " << e.iHH_1 << std::endl;
    std::cout << "iHS: " << log(e.iHH_0/e.iHH_1) << std::endl;
    std::cout << "MAF: " << (double)e.num/(double)hmap.snpLength() << std::endl;
}

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

#include <iostream>

#include "hapbin.hpp"
#include "hapmap.hpp"
#include "ehh.hpp"
#include "argparse.hpp"
#include "ehhpairfinder.hpp"

int main(int argc, char** argv)
{
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<const char*> hap('d', "hap", "Hap file", false, true, "");
    Argument<const char*> map('m', "map", "Map file", false, true, "");
    Argument<double> cutoff('c', "cutoff", "EHH cutoff value (default: 0.05)", false, false, 0.05);
    Argument<double> minMAF('b', "minmaf", "Minimum allele frequency (default: 0.05)", false, false, 0.05);
    Argument<unsigned long long> scale('s', "scale", "Gap scale parameter in bp, used to scale gaps > scale parameter as in Voight, et al.", false, false, 20000);
    Argument<const char*> locus1('1', "locus1", "Locus 1", false, true, 0);
    Argument<const char*> locus2('2', "locus2", "Locus 2", false, true, 0);
    Argument<bool> filternonpoly('f', "filter", "Filter out non-polymorphic loci", false, false);
    Argument<std::string> pops('p', "pop", "Filter by population", true, false,"");
    Argument<std::string> key('k', "key", "Population key", false, false, "");
    ArgParse argparse({&help, &hap, &map, &locus1, &locus2, &cutoff, &minMAF, &scale, &filternonpoly, &pops, &key}, "Usage: ehh2bin --map input.map --hap input.hap --locus1 id1 --locus2 id2");
    argparse.parseArguments(argc, argv);
    HapMap hmap;
    PopKey *popkey = NULL;
    if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
    {
        std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
        return 3;
    }
    if (pops.wasFound() && key.wasFound())
        popkey = new PopKey(key.value(), pops.values());
    if (!hmap.loadHap(hap.value(), filternonpoly.value(), popkey))
    {
        delete popkey;
        return 1;
    }
    hmap.loadMap(map.value());
    EHHPair e;
    std::size_t l1 = hmap.idToIndex(locus1.value());
    std::size_t l2 = hmap.idToIndex(locus2.value());
    if (l1 == std::numeric_limits<std::size_t>::max())
    {
        std::cerr << "no locus with the id: " << locus1.value() << std::endl;
        return 2;
    } else if (l2 == std::numeric_limits<std::size_t>::max())
    {
        std::cerr << "no locus with the id: " << locus2.value() << std::endl;
        return 2;
    }
    EhhPairFinder finder(&hmap, cutoff.value(), minMAF.value(), (double) scale.value(), 1000, true);
    e = finder.calcEhhPair(l1, l2, false);
    e.printEHH(&hmap);
    std::cout << "iHH:    " << e.ihh[0b00] << " " << e.ihh[0b01] << " " << e.ihh[0b10] << " " << e.ihh[0b11] << std::endl;
    std::cout << "iHHNot: " << e.ihhNot[0b00] << " " << e.ihhNot[0b01] << " " << e.ihhNot[0b10] << " " << e.ihhNot[0b11] << std::endl;
    std::cout << "iHS: " << log(e.ihh[0b00]/e.ihhNot[0b00]) << " " << log(e.ihh[0b01]/e.ihhNot[0b01]) << " " << log(e.ihh[0b10]/e.ihhNot[0b10]) << " " << log(e.ihh[0b11]/e.ihhNot[0b11]) << std::endl;
    std::cout << "AF: " << e.af[0b00] << " " << e.af[0b01] << " " << e.af[0b10] << " " << e.af[0b11] << std::endl;
    delete popkey;
}

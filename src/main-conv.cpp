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
#include "argparse.hpp"
#include "ehh.hpp"
#include "ihspairjob.hpp"

#include <functional>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

int tobin(const std::string& in, const std::string& map, const std::string& out, bool filternonpoly, const PopKey* pk, bool filternondiv)
{
    HapMap hap;
    const char* mapname = nullptr;
    if (filternondiv)
        mapname = map.c_str();
    if (!hap.loadHapAscii(in.c_str(), 0, filternonpoly, pk, filternondiv, mapname))
    {
        std::cerr << "Input file " << in << " not found." << std::endl;
        return 1;
    }
    hap.save((out + ".hapbin").c_str());
    if (filternonpoly && map.size() > 0)
        hap.filterMap(map.c_str(), (out + ".map").c_str());
    else if (filternonpoly && map.size() == 0)
        std::cerr << "Warning: Use --map to also filter map." << std::endl;
    std::cout << "Converted haplotype map of dimensions " << hap.snpLength() << "x" << hap.numSnps() << " to binary format." << std::endl;
    return 0;
}

int unstdToTxt(const std::string& in, const std::string& out, const std::string& map, const std::string& hap, bool filter, const PopKey* pk)
{
    HapMap hm;
    hm.loadHap(hap.c_str(),filter,pk);
    hm.loadMap(map.c_str());
    IhsPairJob j;
    j.load(in);
    j.saveAscii(out,hm);
    return 0;
}

int stdToTxt(const std::string& in, const std::string& out, const std::string& map, const std::string& hap, bool filter, const PopKey* pk)
{
    std::ifstream ins(in, std::ios::in | std::ios::binary);
    std::ofstream outs(out);
    uint64_t mnum, count;
    HapMap hm;
    hm.loadHap(hap.c_str(),filter,pk);
    hm.loadMap(map.c_str());
    ins.read((char*) &mnum, sizeof(uint64_t));
    if (mnum != IHS2_STD_MNUM)
    {
        std::cerr << "ERROR: Wrong Magic Number for file " << in << std::endl;
        return 1;
    }
    ins.read((char*) &count, sizeof(uint64_t));
    outs << "focus0 focus1 AF00 StdIHS00 StdNSl00 AF01 StdIHS01 StdNSl01 AF10 StdIHS10 StdNSl10 AF11 StdIHS11 StdNSl11" << std::endl;
    for (uint64_t i = 0; i < count; ++i)
    {
        uint64_t focus0, focus1;
        ins.read((char*) &focus0, sizeof(uint64_t));
        ins.read((char*) &focus1, sizeof(uint64_t));
        outs << hm.indexToId(focus0) << " " << hm.indexToId(focus1) << " ";
        for (std::size_t core = 0b00; core <= 0b11; ++core)
        {
            double af, stdIhs, stdNSl;
            ins.read((char*) &af, sizeof(double));
            ins.read((char*) &stdIhs, sizeof(double));
            ins.read((char*) &stdNSl, sizeof(double));
            outs << af << " " << stdIhs << " " << stdNSl << " ";
        }
        outs << std::endl;
    }
    ins.close();
    outs.close();
    return 0;
}

int main(int argc, char** argv)
{
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<bool> version('v', "version", "Version information", true, false);
    Argument<std::string> hap('d', "hap", "ASCII Hap file", false, false, "");
    Argument<std::string> unstd('u', "unstd", "Unstandardized binary ihs2 file", false, false, "");
    Argument<std::string> std('s', "std", "Standardized binary ihs2 file", false, false, "");
    Argument<bool> filternonpoly('f', "filter-nonpoly", "Filter out non-polymorphic loci", false, false);
    Argument<bool> filternondiv('n', "filter-nondivergent", "Filter out non-divergent loci", false, false);
    Argument<std::string> pops('p', "pop", "Filter by population", true, false,"");
    Argument<std::string> key('k', "key", "Population key", false, false, "");
    Argument<std::string> map('m', "map", "Map file for use with --filter", false, false, "");
    Argument<std::string> outfile('o', "out", "Output file", false, false, "");
    ArgParse argparse({&help, &version, &hap, &unstd, &std, &filternonpoly, &filternondiv, &pops, &key, &map, &outfile}, "Usage: hapbinconv --hap input.hap --out outfile.hapbin\n       hapbinconv --unstd input.dat --out out.txt --hap input.hap --map input.map\n       hapbinconv --std input.std --out out.txt --hap input.hap --map input.map");
    if (!argparse.parseArguments(argc, argv))
        return 1;
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

    std::string out;
    if (outfile.wasFound())
        out = outfile.value();
    if (hap.wasFound() && !(unstd.wasFound() || std.wasFound()))
    {
        if (!outfile.wasFound())
            out = hap.value() + "bin";
        PopKey *popkey = nullptr;
        if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
        {
            std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
            return 3;
        }
        if (pops.wasFound() && key.wasFound())
            popkey = new PopKey(key.value(), pops.values());
        return tobin(hap.value(), map.value(), out, filternonpoly.value(), popkey, filternondiv.value());
    }
    else if (unstd.wasFound())
    {
        if (!map.wasFound() || !hap.wasFound())
        {
            std::cerr << "--unstd requires --hap and --map." << std::endl;
            return 3;
        }
        PopKey *popkey = nullptr;
        if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
        {
            std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
            return 3;
        }
        if (pops.wasFound() && key.wasFound())
            popkey = new PopKey(key.value(), pops.values());
        if (!outfile.wasFound())
            out = unstd.value() + ".txt";
        return unstdToTxt(unstd.value(), out, map.value(), hap.value(), filternonpoly.value(), popkey);
    }
    else if (std.wasFound())
    {
        if (!map.wasFound() || !hap.wasFound())
        {
            std::cerr << "--unstd requires --hap and --map." << std::endl;
            return 3;
        }
        PopKey *popkey = nullptr;
        if ((pops.wasFound() && !key.wasFound()) || (!pops.wasFound() && key.wasFound()))
        {
            std::cerr << "--pop and --key must both be specified to filter by population." << std::endl;
            return 3;
        }
        if (pops.wasFound() && key.wasFound())
            popkey = new PopKey(key.value(), pops.values());
        if (!outfile.wasFound())
            out = std.value() + ".txt";
        return stdToTxt(std.value(), out, map.value(), hap.value(), filternonpoly.value(), popkey);
    }
    return 3;
}


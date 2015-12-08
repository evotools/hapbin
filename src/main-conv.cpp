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
#include "argparse.hpp"

#include <functional>
#include <cstdlib>
#include <iostream>

void tobin(const char* in, const char* out)
{
    HapMap map;
    if (!map.loadHapAscii(in))
    {
        std::cerr << "Input file " << in << " not found." << std::endl;
        return;
    }
    map.save(out);
    std::cout << "Converted haplotype map with " << map.snpLength() << " haplotypes to binary format." << std::endl;
}

int main(int argc, char** argv)
{
    Argument<bool> help('h', "help", "Show this help", true, false);
    Argument<bool> version('v', "version", "Version information", true, false);
    Argument<std::string> hap('d', "hap", "ASCII Hap file", false, false, "");
    Argument<std::string> outfile('o', "out", "Binary output file", false, false, "out.hapbin");
    ArgParse argparse({&help, &version, &hap, &outfile}, "Usage: hapbinconv --hap input.hap --out outfile.hapbin");
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
    else if (!hap.wasFound()) {
        std::cout << "Please specify --hap." << std::endl;
        return 2;
    }
    tobin(hap.value().c_str(), outfile.value().c_str()); 
    return 0;
}



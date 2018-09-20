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

#include "config.h"
#include <string>

void calcIhsNoMpi(
    const std::string& hap,
    const std::string& map,
    const std::string& outfile,
    double cutoff,
    double minMAF,
    double scale,
    unsigned long long maxExtend,
    int bins,
    bool binom);

void calcIhsMpi(
    const std::string& hapfile,
    const std::string& mapfile,
    const std::string& outfile,
    double cutoff,
    double minMAF,
    double scale,
    unsigned long long maxExtend,
    int binFactor,
    bool binom);

void calcXpehhNoMpi(
    const std::string& hapA,
    const std::string& hapB,
    const std::string& map,
    const std::string& outfile,
    double cutoff,
    double minMAF,
    double scale,
    unsigned long long maxExtend,
    int bins,
    bool binom);

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
    bool binom);

#if MPI_FOUND
class ParameterStream;
struct IhsScore;
struct XPEHH;
ParameterStream& operator<<(ParameterStream& out, const IhsScore& info);
ParameterStream& operator>>(ParameterStream& in, IhsScore& info);
ParameterStream& operator<<(ParameterStream& out, const XPEHH& info);
ParameterStream& operator>>(ParameterStream& in, XPEHH& info);

#define calcIhs calcIhsMpi
#define calcXpehh calcXpehhMpi
#else
#define calcIhs calcIhsNoMpi
#define calcXpehh calcXpehhNoMpi
#endif


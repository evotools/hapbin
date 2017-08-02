/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2015-2017 Colin MacLean <cmaclean@illinois.edu>
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

#include "ihspairjob.hpp"
#include "hapmap.hpp"
#include <fstream>
#include <cmath>

IhsPairJob::IhsPairJob(std::size_t start, std::size_t end)
    : m_range{start,end}
    , m_sumIhs{}
    , m_results{}
{}

IhsPairJob::IhsPairJob()
    : m_range{}
    , m_sumIhs{}
    , m_results{}
{}

void IhsPairJob::add(const EHHPair& result)
{
    for (std::size_t i = 0b00; i <= 0b11; ++i)
    {
        double ihs = log(result.ihh[i]/result.ihhNot[i]);
        if (!std::isnan(ihs) && !std::isinf(ihs))
            m_sumIhs[i] += ihs;
    }
    m_results.push_back(std::move(result));
}

void IhsPairJob::save(const std::string& file)
{
    std::ofstream f(file, std::ios::out | std::ios::binary);
    if (!f.good())
    {
        std::cerr << "ERROR: Could not open " << file << " for writing!" << std::endl;
    }
    f.write((char*) &IhsPairJob::magicNumber, sizeof(uint64_t));
    uint64_t start = m_range[0];
    uint64_t end   = m_range[1];
    f.write((char*) &start, sizeof(uint64_t));
    f.write((char*) &end, sizeof(uint64_t));
    uint64_t num = m_results.size();
    f.write((char*) &num, sizeof(uint64_t));
    f.write((char*) &m_sumIhs, sizeof(double));
    for (const auto& it : m_results)
    {
        uint64_t focus0 = it.focus[0];
        uint64_t focus1 = it.focus[1];
        f.write((char*) &focus0, sizeof(uint64_t));
        f.write((char*) &focus1, sizeof(uint64_t));
        for(std::size_t core = 0b00; core <= 0b11; ++core)
        {
            double af = it.af[core];
            double ihh = it.ihh[core];
            double ihhNot = it.ihhNot[core];
            double sl = it.sl[core];
            double slNot = it.slNot[core];
            f.write((char*) &af, sizeof(decltype(af)));
            f.write((char*) &ihh, sizeof(decltype(ihh)));
            f.write((char*) &ihhNot, sizeof(decltype(ihhNot)));
            f.write((char*) &sl, sizeof(decltype(sl)));
            f.write((char*) &slNot, sizeof(decltype(slNot)));
        }
    }
}

void IhsPairJob::saveAscii(const std::string& out, const HapMap& hm)
{
    std::ofstream f(out);
    if (!f.good())
    {
        std::cerr << "ERROR: Could not open " << out << " for writing!" << std::endl;
    }
    f << "focus0 focus1 AF00\tIHH00\tIHHNot00\tIHS00\tSL00\tSLNot00\tnSL00 "
                       "AF01\tIHH01\tIHHNot01\tIHS01\tSL01\tSLNot01\tnSL01 "
                       "AF10\tIHH10\tIHHNot10\tIHS10\tSL10\tSLNot10\tnSL10 "
                       "AF11\tIHH11\tIHHNot11\tIHS11\tSL11\tSLNot11\tnSL11" << std::endl;;
    for (const auto& it : m_results)
    {
        f << hm.indexToId(it.focus[0]) << "\t" << hm.indexToId(it.focus[1]) << "\t";
        for(std::size_t core = 0b00; core <= 0b11; ++core)
        {
            f << it.af[core] << "\t";
            double ihh = it.ihh[core];
            double ihhNot = it.ihhNot[core];
            double sl = it.sl[core];
            double slNot = it.slNot[core];
            f << ihh << "\t" << ihhNot << "\t" << log(ihh/ihhNot) << "\t";
            f << sl  << "\t" << slNot  << "\t" << log( sl/slNot );
            if (core == 0b11)
                f << std::endl;
            else
                f << "\t";
        }
    }
}

int IhsPairJob::load(const std::string& file)
{
    std::ifstream f(file, std::ios::in | std::ios::binary);
    if (!f.good())
    {
        //std::cerr << "ERROR: Could not open result file " << file << "!" << std::endl;
        return 1;
    }
    uint64_t mN;
    f.read((char*) &mN, sizeof(uint64_t));
    if (mN != magicNumber)
    {
        std::cerr << "ERROR: Bad magic number for file " << file << "!" << std::endl;
        return 2;
    }
    uint64_t start, end, num;
    try {
        f.exceptions(std::ifstream::eofbit);
        f.read((char*) &start, sizeof(uint64_t));
        m_range[0] = start;
        f.read((char*) &end, sizeof(uint64_t));
        m_range[1] = end;
        f.read((char*) &num, sizeof(uint64_t));
        f.read((char*) &m_sumIhs, sizeof(double));
        for (std::size_t i = 0; i < num; ++i)
        {
            EHHPair p;
            uint64_t focus0, focus1;
            f.read((char*) &focus0, sizeof(decltype(focus0)));
            f.read((char*) &focus1, sizeof(decltype(focus1)));
            p.focus[0] = focus0;
            p.focus[1] = focus1;
            for(std::size_t core = 0b00; core <= 0b11; ++core)
            {
                f.read((char*) &p.af[core], sizeof(decltype(p.af[0])));
                f.read((char*) &p.ihh[core], sizeof(decltype(p.ihh[0])));
                f.read((char*) &p.ihhNot[core], sizeof(decltype(p.ihhNot[0])));
                f.read((char*) &p.sl[core], sizeof(decltype(p.sl[0])));
                f.read((char*) &p.slNot[core], sizeof(decltype(p.slNot[0])));
            }
            m_results.push_back(std::move(p));
        }
        return 0;
    }
    catch (std::ios_base::failure &e)
    {
        std::cout << "Recalculating " << start << "-" << end << ". Truncated file." << std::endl;
        return 3;
    }
}

std::vector<EHHPair> IhsPairJob::results() const
{
    return m_results;
}

const uint64_t IhsPairJob::magicNumber = IHS2_JOB_MNUM;

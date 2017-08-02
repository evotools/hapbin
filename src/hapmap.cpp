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
#include <cassert>
#include <fstream>
#include <cstdlib>
#include <algorithm>

HapMap::HapMap()
    : m_idMap{}
    , m_physPos(nullptr)
    , m_genPos(nullptr)
    , m_numSnps{}
    , m_data(nullptr)
    , m_snpLength{}
    , m_snpDataSize{}
    , m_snpDataSize64{}
    , m_snpDataSizeULL{}
    , m_polySites{NULL}
    , m_polyCount{0}
{

}

HapMap::~HapMap()
{
    aligned_free(m_data);
    delete[] m_physPos;
    delete[] m_genPos;
}

void HapMap::loadMap(const char* mapFilename)
{
    assert(m_numSnps != 0ULL); //Must loadHap first.
    if (m_physPos)
        delete[] m_physPos;
    if (m_genPos)
        delete[] m_genPos;
    m_physPos = new unsigned long long[m_numSnps];
    m_genPos = new double[m_numSnps];

    std::ifstream file(mapFilename);

    std::string line;
    std::size_t lineIndex = 0ULL, mapIndex = 0ULL;

    while (std::getline(file, line))
    {
        if (mapIndex == m_numSnps)
        {
            std::cerr << "WARNING: Map file has more loci than hap file!" << std::endl;
            return;
        }
        std::vector<std::string> split = splitString(line, ' ');
        if (split.size() == 1)
            split = splitString(line, '\t');
        if (split.size() == 4)
        {
            m_idMap[mapIndex] = split[1];
            m_genPos[mapIndex] = atof(split[2].c_str());
            m_physPos[mapIndex] = strtoull(split[3].c_str(), 0, 10);
            if (m_polyCount == m_numSnps)
            {
                if (m_polySites[lineIndex])
                    ++mapIndex;
            }
            else
            {
                ++mapIndex;
            }
            ++lineIndex;

        }
    }
    if (mapIndex != m_numSnps)
    {
        std::cerr<< "ERROR: Map file must have the same number of loci as hap file! Hap file has " << m_numSnps << " SNPs. Map file has " << mapIndex << " SNPs." << std::endl;
        abort();
    }
    m_polySites.resize(0);
    file.close();
}

void HapMap::filterMap(const char* inmap, const char* outmap)
{
    assert(m_polyCount == m_numSnps && m_polySites.size() > 0);

    std::ifstream in(inmap);
    std::ofstream out(outmap);

    std::string line;

    std::size_t lineNum = 0ULL;
    while (std::getline(in, line))
    {
        if (m_polySites[lineNum])
            out << line << std::endl;
        ++lineNum;
    }
}

void HapMap::save(const char* filename)
{
    std::ofstream f(filename, std::ios::out | std::ios::binary);
    f.write((char*) &magicNumber, sizeof(uint64_t));
    uint64_t ns = m_numSnps;
    f.write((char*) &ns, sizeof(uint64_t));
    uint64_t sl = this->snpLength();
    f.write((char*) &sl, sizeof(uint64_t));

    for(uint64_t i = 0; i < m_numSnps; ++i)
        f.write((char*) &m_data[i*m_snpDataSize], m_snpDataSize64*sizeof(uint64_t));
    f.close();
}

void HapMap::setSnpLength(uint64_t l)
{
    m_snpLength = l;
    m_snpDataSize = ::bitsetSize<PrimitiveType>(l);
    m_snpDataSizeULL = ::bitsetSize<unsigned long long>(l);
    m_snpDataSize64 = ::bitsetSize<uint64_t>(l);
}

unsigned long long HapMap::physicalPosition(std::size_t index) const
{
    return m_physPos[index];
}

double HapMap::geneticPosition(std::size_t index) const
{
    return m_genPos[index];
}

std::string HapMap::indexToId(std::size_t index) const
{
    return m_idMap.at(index);
}

std::size_t HapMap::idToIndex(const std::string& id) const
{
    for(const auto &p : m_idMap)
        if(p.second == id)
            return p.first;
    return std::numeric_limits<std::size_t>::max();
}

std::size_t HapMap::querySnpLength(const char* filename)
{
    std::ifstream f(filename, std::ios::in | std::ios::binary);
    if (!f.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }
    uint64_t check;
    f.read((char*) &check, sizeof(uint64_t));
    f.close();
    if (check == magicNumber)
        return querySnpLengthBinary(filename);
    else
        return querySnpLengthAscii(filename);
}

std::size_t HapMap::querySnpLengthBinary(const char* filename)
{
    std::ifstream f(filename, std::ios::in | std::ios::binary);

    uint64_t sl;
    f.seekg(sizeof(uint64_t)*2);
    f.read((char*) &sl, sizeof(uint64_t));
    return sl;
}

std::size_t HapMap::querySnpLengthAscii(const char* filename)
{
    std::ifstream f(filename);
    std::string line;
    std::getline(f, line);
    std::size_t sl = (line.size()+1)/2;
    return sl;
}

bool HapMap::loadHap(const char* filename, bool onlyPoly, const PopKey* popkey)
{
    std::ifstream f(filename, std::ios::in | std::ios::binary);
    if (!f.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }
    uint64_t check;
    f.read((char*) &check, sizeof(uint64_t));
    f.close();
    if (check == magicNumber)
        return loadHapBinary(filename, onlyPoly, popkey);
    else
        return loadHapAscii(filename, 0, onlyPoly, popkey);
}

bool HapMap::loadHapBinary(const char* filename, bool onlyPoly, const PopKey* popkey, bool removeNonDiverging, const char* mapFilename)
{
    aligned_free(m_data);

    std::ifstream f(filename, std::ios::in | std::ios::binary);
    std::ifstream *mapfile = nullptr;
    if (mapFilename)
        mapfile = new std::ifstream(mapFilename);
    if (!f.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }

    uint64_t check;
    f.read((char*) &check, sizeof(uint64_t));
    if (check != magicNumber)
    {
        std::cerr << "ERROR: Wrong file type: " << filename << ". Expected binary format." << std::endl;
        return false;
    }

    f.read((char*) &m_numSnps, sizeof(uint64_t));
    f.read((char*) &m_snpLength, sizeof(uint64_t));

    m_snpDataSize = ::bitsetSize<PrimitiveType>(m_snpLength);
    m_snpDataSize64 = ::bitsetSize<uint64_t>(m_snpLength);
    m_snpDataSizeULL = ::bitsetSize<unsigned long long>(m_snpLength);

    std::size_t origDataSize, origDataSize64, origDataSizeULL, origLength;
    if (popkey != NULL)
    {
        origDataSize = m_snpDataSize;
        origDataSize64 = m_snpDataSize64;
        origDataSizeULL = m_snpDataSizeULL;
        origLength = m_snpLength;
        m_snpLength = popkey->count();
        m_snpDataSize = ::bitsetSize<PrimitiveType>(m_snpLength);
        m_snpDataSize64 = ::bitsetSize<uint64_t>(m_snpLength);
        m_snpDataSizeULL = ::bitsetSize<unsigned long long>(m_snpLength);
    }
    m_data = (PrimitiveType*) aligned_alloc(128, m_snpDataSize*m_numSnps*sizeof(PrimitiveType));

    m_polySites.resize(0);
    m_polySites.reserve(m_numSnps);
    m_polyCount = 0ULL;
    std::size_t index = 0ULL;
    if (onlyPoly && popkey != NULL)
    {
        unsigned long long* buffer = new unsigned long long[origDataSize64]{};

        double genPos = -1.0;
        double prevGenPos = -1.0;
        for(uint64_t i = 0; i < m_numSnps; ++i)
        {
            unsigned long long* snpULL = (unsigned long long*) &m_data[m_snpDataSize*index];
            unsigned long long* prevSnpULL = (unsigned long long*) &m_data[m_snpDataSize*(index-1)];
            std::fill_n(snpULL, m_snpDataSize*VEC, 0ULL);
            f.read((char*) buffer, sizeof(uint64_t)*origDataSize64);
            if (removeNonDiverging && mapfile)
            {
                do
                {
                    std::string line;
                    if (!std::getline(*mapfile, line))
                    {
                        std::cerr << "Bad map file length!" << std::endl;
                        break;
                    }
                    std::vector<std::string> split = splitString(line, ' ');
                    if (split.size() == 1)
                        split = splitString(line, '\t');
                    if (split.size() == 4)
                    {
                        genPos = atof(split[2].c_str());
                    }
                    else
                    {
                        continue;
                    }
                } while (0);
            }
            filter(buffer, snpULL, origLength, m_snpLength, *popkey);
            bool diverges = false;
            std::size_t count = 0;
            for(std::size_t j = 0; j < m_snpDataSizeULL; ++j)
            {
                if (index > 0 && removeNonDiverging)
                    diverges |= (bool) (snpULL[j] ^ prevSnpULL[j]);
                count += popcount1(snpULL[j]);
            }
            diverges = (diverges && (genPos != prevGenPos)) || !removeNonDiverging;
            prevGenPos = genPos;
            if (count == 0 || count == m_snpLength || !diverges)
            {
                m_polySites.push_back(false);
            }
            else
            {
                m_polySites.push_back(true);
                ++m_polyCount;
                ++index;
            }
        }
        delete[] buffer;
        m_numSnps = index;
    }
    else if (onlyPoly && popkey == NULL)
    {
        double genPos = -1.0;
        double prevGenPos = -1.0;
        for(uint64_t i = 0; i < m_numSnps; ++i)
        {
            uint64_t* snp64 = (uint64_t*) &m_data[m_snpDataSize*index];
            unsigned long long* snpULL = (unsigned long long*) snp64;
            unsigned long long* prevSnpULL = (unsigned long long*) &m_data[m_snpDataSize*(index-1)];
            f.read((char*) &m_data[i*this->m_snpDataSize], sizeof(uint64_t)*m_snpDataSize64);
            if (removeNonDiverging && mapfile)
            {
                do
                {
                    std::string line;
                    if (!std::getline(*mapfile, line))
                    {
                        std::cerr << "Bad map file length!" << std::endl;
                        break;
                    }
                    std::vector<std::string> split = splitString(line, ' ');
                    if (split.size() == 1)
                        split = splitString(line, '\t');
                    if (split.size() == 4)
                    {
                        genPos = atof(split[2].c_str());
                    }
                    else
                    {
                        continue;
                    }
                } while (0);
            }
            bool diverges = false;
            std::size_t count = 0;
            for(std::size_t j = 0; j < m_snpDataSizeULL; ++j)
            {
                if (index > 0 && removeNonDiverging)
                    diverges |= (bool) (snpULL[j] ^ prevSnpULL[j]);
                count += popcount1(snpULL[j]);
            }
            diverges = (diverges && (genPos != prevGenPos)) || !removeNonDiverging;
            prevGenPos = genPos;
            if (count == 0 || count == m_snpLength || !diverges)
            {
                m_polySites.push_back(false);
            }
            else
            {
                m_polySites.push_back(true);
                ++m_polyCount;
                ++index;
            }
        }
        m_numSnps = index;
    }
    else if (!onlyPoly && popkey != NULL)
    {
        unsigned long long* buffer = new unsigned long long[origDataSizeULL]{};
        for(uint64_t i = 0; i < m_numSnps; ++i)
        {
            f.read((char*) buffer, sizeof(uint64_t)*origDataSize64);
            filter(buffer, (unsigned long long*) &m_data[m_snpDataSize*i], origLength, m_snpLength, *popkey);
        }
        delete[] buffer;
    }
    else
    {
        for(uint64_t i = 0; i < m_numSnps; ++i)
        {
            f.read((char*) &m_data[i*this->m_snpDataSize], sizeof(uint64_t)*m_snpDataSize64);
            ++index;
        }
    }

    return true;
}

void HapMap::filterNonDiverging(bool keepGenPosDivergence)
{
    for (std::size_t i = 1; i < m_numSnps; ++i)
    {
        unsigned long long* snpULL = (unsigned long long*) &m_data[m_snpDataSize*i];
        unsigned long long* prevSnpULL = (unsigned long long*) &m_data[m_snpDataSize*(i-1)];
        bool diverges = false;
        for (std::size_t j = 0; j < m_snpDataSizeULL; ++j)
        {
            diverges |= (bool) (snpULL[j] ^ prevSnpULL[j]);
        }
        if (keepGenPosDivergence)
            diverges = diverges && (m_genPos[i] == m_genPos[i-1]);
    }
}

bool HapMap::loadHapAscii(const char* filename, std::size_t maxLength, bool onlyPoly, const PopKey* popkey, bool removeNonDiverging, const char* mapFilename)
{
    aligned_free(m_data);

    std::ifstream file(filename);
    std::ifstream *mapfile = nullptr;
    if (mapFilename)
        mapfile = new std::ifstream(mapFilename);
    if (!file.good())
    {
        std::cerr << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }
    m_numSnps = 1ULL;
    std::string line;

    std::getline(file, line);
    m_snpLength = (line.size()+1)/2;
    std::size_t origSnpLength = m_snpLength;
    if (maxLength  > 0 && maxLength < origSnpLength)
        origSnpLength = maxLength;
    if(popkey)
        m_snpLength = popkey->count();
    if (maxLength  > 0 && maxLength < m_snpLength)
        m_snpLength = maxLength;

    //while(std::getline(file, line))
    //    ++m_numSnps;

    m_snpDataSize = ::bitsetSize<PrimitiveType>(m_snpLength);
    m_snpDataSize64 = ::bitsetSize<uint64_t>(m_snpLength);
    m_snpDataSizeULL = ::bitsetSize<unsigned long long>(m_snpLength);

    std::vector<unsigned long long> vec;

    file.clear();
    file.seekg(0, std::ios::beg);

    if (onlyPoly)
    {
        m_polySites.resize(0);
    }

    std::size_t index = 0;
    double genPos = -1.0;
    double prevGenPos = -1.0;
    std::size_t numNonDiverging = 0;
    unsigned long long* prevSnp = new unsigned long long[m_snpDataSize*VEC]{};
    //std::cout << "have it: " << removeNonDiverging << " " << mapfile << std::endl;
    std::cout << "SNPs read: 0";
    for (std::size_t i = 0; std::getline(file, line); ++i)
    {

        if (line.size() > 0)
        {
            vec.resize((index+1)*m_snpDataSize*VEC, 0);
            unsigned long long* snp = &vec[m_snpDataSize*VEC*index];
            if (popkey)
                convert<unsigned long long>(line.c_str(), snp, *popkey, maxLength);
            else
                convert<unsigned long long>(line.c_str(), snp, maxLength);
            if(onlyPoly)
            {
                if (removeNonDiverging && mapfile)
                {
                    do
                    {
                        std::string line;
                        if (!std::getline(*mapfile, line))
                        {
                            std::cerr << "Bad map file length!" << std::endl;
                            break;
                        }
                        std::vector<std::string> split = splitString(line, ' ');
                        if (split.size() == 1)
                            split = splitString(line, '\t');
                        if (split.size() == 4)
                        {
                            genPos = atof(split[2].c_str());
                        }
                        else
                        {
                            continue;
                        }
                    } while (0);
                }
                bool diverges = false;
                unsigned long long divergence_tracker = 0ULL;
                std::size_t count = 0;
                for(std::size_t j = 0; j < m_snpDataSizeULL; ++j)
                {
                    if (index > 0 && removeNonDiverging)
                    {
                        //std::cout << snp[j] << " " << prevSnp[j] << std::endl;
                        divergence_tracker |= (snp[j] ^ prevSnp[j]);
                    }
                    count += popcount1(snp[j]);
                }
                //std::cout << std::endl;
                diverges = divergence_tracker;
                if (index == 0)
                    diverges = true;
                //print_snp(snp,m_snpLength);
                //std::cout << "diverges: " << divergence_tracker << " " << diverges << " " << count << std::endl;
                diverges = (diverges && (genPos != prevGenPos)) || !removeNonDiverging;
                std::copy_n(snp,m_snpDataSizeULL,prevSnp);
                prevGenPos = genPos;
                if (!diverges && count > 0 && count < m_snpLength)
                    ++numNonDiverging;
                if (count == 0 || count == m_snpLength || !diverges)
                {
                    m_polySites.push_back(false);
                    std::fill_n(snp, m_snpDataSize*VEC, 0LL);
                }
                else
                {
                    m_polySites.push_back(true);
                    ++m_polyCount;
                    ++index;
                }
            }
            else
            {
                ++index;
            }
            if (index % 1000 == 0)
                std::cout << "\rSNPs read: " << index;
        }
    }
    m_numSnps = index;
    std::cout << "\rSNPs read: " << index << std::endl;
    std::cout << "Number of non-diverging sites: " << numNonDiverging << std::endl;
    /*
     * Allocate memory aligned to 128 bytes (cache line). Must be aligned to at least 32 bytes for required AVX instructions.
     */
    m_data = (PrimitiveType*) aligned_alloc(128, m_snpDataSize*m_numSnps*sizeof(PrimitiveType));
    std::copy(vec.cbegin(), vec.cend(), (unsigned long long*) m_data);
    return true;
}

std::size_t HapMap::count(std::size_t focus1) const
{
    int count = 0;
    unsigned long long* hapsULL = (unsigned long long*) &m_data[focus1*m_snpDataSize];
    for(std::size_t i = 0; i < m_snpDataSizeULL; ++i)
    {
        count += popcount1(hapsULL[i]);
    }
    return count;
}

const HapMap::PrimitiveType* HapMap::operator[](std::size_t i) const
{
    return &m_data[i*m_snpDataSize];
}

const uint64_t HapMap::magicNumber = 3544454305642733928ULL;

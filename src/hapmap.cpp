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

#include "hapmap.hpp"
#include <cassert>
#include <cstdlib>

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
{

}

HapMap::~HapMap()
{
    aligned_free(m_data);
    delete m_physPos;
    delete m_genPos;
}

void HapMap::loadMap(const char* mapFilename)
{
    assert(m_numSnps != 0ULL); //Must loadHap first.
    if (m_physPos)
        delete m_physPos;
    if (m_genPos)
        delete m_genPos;
    m_physPos = new unsigned long long[m_numSnps];
    m_genPos = new double[m_numSnps];
        
    std::ifstream file(mapFilename);

    std::string line;
    
    std::size_t lineNum = 0ULL;
    while (std::getline(file, line)) 
    {
        if (lineNum == m_numSnps)
        {
            std::cerr << "WARNING: Map file has more loci than hap file!" << std::endl;
            return;
        }
        std::vector<std::string> split = splitString(line, ' ');
	if (split.size() == 1)
            split = splitString(line, '\t');
        if (split.size() == 4) {
            m_idMap[lineNum] = split[1];
            m_genPos[lineNum] = atof(split[2].c_str());
            m_physPos[lineNum] = strtoull(split[3].c_str(), 0, 10);
            ++lineNum;
        }
    }
    if (lineNum != m_numSnps)
    {
        std::cerr << "ERROR: Map file must have the same number of loci as hap file! Hap file has " << m_numSnps << " SNPs. Map file has " << lineNum << " SNPs." << std::endl;
        if (lineNum == 0)
            std::cerr << "Perhaps the Map file format is wrong?" << std::endl;
        abort();
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
}

void HapMap::setSnpLength(uint64_t l)
{
    m_snpLength = l; 
    m_snpDataSize = ::bitsetSize<PrimitiveType>(l);
    m_snpDataSizeULL = ::bitsetSize<unsigned long long>(l);
    m_snpDataSize64 = ::bitsetSize<uint64_t>(l);
}

unsigned long long HapMap::physicalPosition(std::size_t line) const
{
    return m_physPos[line];
}

double HapMap::geneticPosition(std::size_t line) const
{
    return m_genPos[line];
}

std::string HapMap::lineToId(std::size_t line) const
{
    return m_idMap.at(line);
}

std::size_t HapMap::idToLine(const std::string& id) const
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
    f.seekg(sizeof(uint64_t));
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

bool HapMap::loadHap(const char* filename)
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
        return loadHapBinary(filename);
    else
        return loadHapAscii(filename);
}

bool HapMap::loadHapBinary(const char* filename)
{
    aligned_free(m_data);
    
    std::ifstream f(filename, std::ios::in | std::ios::binary);
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
    m_data = (PrimitiveType*) aligned_alloc(128, m_snpDataSize*m_numSnps*sizeof(PrimitiveType));
    
    for(uint64_t i = 0; i < m_numSnps; ++i)
        f.read((char*) &m_data[i*this->m_snpDataSize], sizeof(uint64_t)*m_snpDataSize64);
    
    return true;
}

bool HapMap::loadHapAscii(const char* filename, std::size_t maxLength)
{
    aligned_free(m_data);
    
    std::ifstream file(filename);
    if (!file.good())
    {
        std::cout << "ERROR: Cannot open file or file not found: " << filename << std::endl;
        return false;
    }
    m_numSnps = 1ULL;
    std::string line;
    
    std::getline(file, line);
    m_snpLength = (line.size()+1)/2;
    if (maxLength  > 0 && maxLength < m_snpLength)
        m_snpLength = maxLength;
    while(std::getline(file, line))
        ++m_numSnps;
    
    m_snpDataSize = ::bitsetSize<PrimitiveType>(m_snpLength);
    m_snpDataSize64 = ::bitsetSize<uint64_t>(m_snpLength);
    m_snpDataSizeULL = ::bitsetSize<unsigned long long>(m_snpLength);
    
    /*
     * Allocate memory aligned to 128 bytes (cache line). Must be aligned to at least 32 bytes for required AVX instructions.
     */
    m_data = (PrimitiveType*) aligned_alloc(128, m_snpDataSize*m_numSnps*sizeof(PrimitiveType));
    
    for (size_t i = 0; i < this->m_snpDataSize*this->m_numSnps; ++i)
    {
        m_data[i] = ZERO;
    }
    
    file.clear();
    file.seekg(0, std::ios::beg);
    
    for (size_t i = 0; i < m_numSnps; ++i)
    {
        std::getline(file, line);
        if (line.size() > 0)
        {
            convert<unsigned long long>(line.c_str(), (unsigned long long*) (&m_data[this->m_snpDataSize*i]), maxLength);
        }
    }
    return true;
}

const uint64_t HapMap::magicNumber = 3544454305642733928ULL;

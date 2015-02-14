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

#include <cstdint>
#include <vector>
#include <string>
#include <bitset>
#include <fstream>
#include <unordered_map>
#include <map>
#include "hapbin.hpp"

#ifndef HAPMAP_HPP
#define HAPMAP_HPP

#include <stdexcept>
#include <iostream>

class HapMap
{
public:
#if VEC==4
    using PrimitiveType = v4ul;
#elif VEC==2
    using PrimitiveType = v2ul;
#else
    using PrimitiveType = unsigned long long;
#endif
    
    HapMap();
    static std::size_t querySnpLengthBinary(const char* filename);
    static std::size_t querySnpLengthAscii(const char* filename);
    static std::size_t querySnpLength(const char* filename);
    double geneticPosition(std::size_t line) const;
    long long unsigned int physicalPosition(std::size_t line) const;
    void loadMap(const char* mapFileName);
    std::string lineToId(std::size_t line) const;
    std::size_t idToLine(const std::string& id) const;
    bool loadHapBinary(const char* filename);
    bool loadHapAscii(const char* filename, std::size_t maxLength = 0);
    bool loadHap(const char* filename);
    void save(const char* filename);
    
    std::size_t numSnps() const { return m_numSnps; }
    std::size_t snpLength() const { return m_snpLength; }
    void setSnpLength(std::size_t l);
    std::size_t snpDataSize() const { return m_snpDataSize; }
    std::size_t snpDataSizeULL() const { return m_snpDataSizeULL; }
    std::size_t snpDataSize64() const { return m_snpDataSize64; }
    PrimitiveType* rawData() { return m_data; }
    ~HapMap();
    
    static const uint64_t magicNumber;
    
protected:
    std::map<std::size_t, std::string> m_idMap;
    unsigned long long* m_physPos;
    double* m_genPos;
    
    std::size_t m_numSnps;
    PrimitiveType *m_data;
    std::size_t m_snpLength;
    std::size_t m_snpDataSize;
    std::size_t m_snpDataSize64;
    std::size_t m_snpDataSizeULL;
};

#endif // CTCHAPM_HPP

/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2016-2017 Colin MacLean <cmaclean@illinois.edu>
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

#include "popkey.hpp"
#include "hapbin.hpp"
#include <fstream>
#include <algorithm>
#include <iostream>

PopKey::PopKey(const std::string& keyfilename, const std::vector<std::string>& pops)
    : m_count{}
{
    std::ifstream keyfile(keyfilename);
    std::string line;

    std::getline(keyfile,line); //skip header

    while (std::getline(keyfile, line)) {
        auto s = splitString(line,' ');
        auto pop = s[1];
        auto group = s[2];
        auto it = std::find(pops.cbegin(), pops.cend(), pop);
        auto it2 = std::find(pops.cbegin(), pops.cend(), group);
        if (it != pops.cend() || it2 != pops.cend())
        {
            m_key.push_back(true);
            m_key.push_back(true);
            m_count+=2;
        }
        else
        {
            m_key.push_back(false);
            m_key.push_back(false);
        }
    }
}

PopKey::PopKey(const char* mask)
    : m_count{}
{
    std::size_t i = 0;
    char k = mask[i];
    while (k != '\0')
    {
        switch (k)
        {
            case '1':
                m_key.push_back(true);
                ++m_count;
                break;
            case '0':
                m_key.push_back(false);
                break;
            case ' ':
                break;
            default:
                std::cerr << "Error reading key: '" << k << "'" << std::endl;
                abort();
        }
        ++i;
        k = mask[i];
    }
}

bool PopKey::operator[](std::size_t index) const
{
    return m_key.at(index);
}

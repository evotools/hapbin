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

#ifndef POPKEY_H
#define POPKEY_H

#include <string>
#include <vector>
#include <cstdint>

class PopKey
{
public:
    PopKey(const std::string& keyfilename, const std::vector<std::string>& pops);
    PopKey(const char* mask);
    bool operator[](std::size_t index) const;
    std::size_t count() const { return m_count; }
protected:
    std::vector<bool> m_key;
    std::size_t m_count;
};

#endif // POPKEY_H

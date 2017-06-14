/*
 * MPIRPC: MPI based invocation of functions on other ranks
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

#ifndef PARAMETERSTREAM_H
#define PARAMETERSTREAM_H

#include<vector>
#include<cstddef>
#include<cstdint>
#include<string>
#include<sstream>
#include<map>

class ParameterStream
{
public:
    ParameterStream() = delete;
    ParameterStream(std::vector<char>* buffer);
    //ParameterStream(const char* data, size_t length);

    void seek(std::size_t pos);

    void writeBytes(const char* b, size_t length);
    void readBytes(char*& b, size_t length);
    char* data();
    const char* constData() const;
    std::vector<char>* dataVector() const;
    size_t size() const;

    ParameterStream& operator<<(int8_t val);
    ParameterStream& operator<<(int16_t val);
    ParameterStream& operator<<(int32_t val);
    ParameterStream& operator<<(int64_t val);
    ParameterStream& operator<<(uint8_t val);
    ParameterStream& operator<<(uint16_t val);
    ParameterStream& operator<<(uint32_t val);
    ParameterStream& operator<<(uint64_t val);
    ParameterStream& operator<<(long long val);
    ParameterStream& operator<<(unsigned long long val);

    ParameterStream& operator>>(int8_t& val);
    ParameterStream& operator>>(int16_t& val);
    ParameterStream& operator>>(int32_t& val);
    ParameterStream& operator>>(int64_t& val);
    ParameterStream& operator>>(uint8_t& val);
    ParameterStream& operator>>(uint16_t& val);
    ParameterStream& operator>>(uint32_t& val);
    ParameterStream& operator>>(uint64_t& val);
    ParameterStream& operator>>(long long& val);
    ParameterStream& operator>>(unsigned long long& val);

    ParameterStream& operator<<(float val);
    ParameterStream& operator<<(double val);
    ParameterStream& operator>>(float& val);
    ParameterStream& operator>>(double& val);

    ParameterStream& operator<<(bool val);
    ParameterStream& operator>>(bool& val);

    ParameterStream& operator<<(const char* s);
    ParameterStream& operator>>(char *& s);

    ParameterStream& operator<<(const std::string& val);
    ParameterStream& operator>>(std::string& val);

protected:
    std::vector<char> *m_data;
    std::size_t  m_pos;
};

template<typename T>
inline void marshal(ParameterStream& s, T val) {
    s << val;
}

template<typename T>
inline typename std::decay<T>::type unmarshal(ParameterStream& s) {
    typename std::decay<T>::type ret;
    s >> ret;
    return ret;
}

template<typename T, typename... O>
ParameterStream& operator<<(ParameterStream& out, const std::vector<T,O...>& vector)
{
    out << vector.size();
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
        out << vector[i];
    }
    return out;
}

template <typename T, typename... O>
ParameterStream& operator>>(ParameterStream& in, std::vector<T, O...>& vector)
{
    std::size_t size;
    in >> size;
    vector.resize(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        T val;
        in >> val;
        vector[i] = val;
    }
    return in;
}

template<typename T, typename U, typename... O>
ParameterStream& operator<<(ParameterStream& out, const std::map<T, U, O...>& map)
{
    out <<  map.size();
    for (auto pair : map)
    {
        out << pair.first << pair.second;
    }
    return out;
}

template<typename T, typename U, typename... O>
ParameterStream& operator>>(ParameterStream& in, std::map<T, U, O...>& map)
{
    std::size_t size;
    in >> size;
    for (std::size_t i = 0; i < size; ++i)
    {
        T first;
        U second;
        in >> first >> second;
        map.insert(std::make_pair(first, second));
    }
    return in;
}

template<typename T, typename U>
ParameterStream& operator<<(ParameterStream& out, const std::pair<T,U>& p)
{
    out << p.first << p.second;
    return out;
}

template<typename T, typename U>
ParameterStream& operator>>(ParameterStream& in, std::pair<T,U>& p)
{
    in >> p.first >> p.second;
    return in;
}

/*template<typename... Ts>
ParameterStream& operator<<(ParameterStream& out, const std::tuple<Ts...>& t)
{
    out <<
}*/

#endif // PARAMETERSTREAM_H

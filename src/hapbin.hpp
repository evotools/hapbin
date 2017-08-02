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

#ifndef HAPBIN_HPP
#define HAPBIN_HPP

#include <type_traits>
#include <cassert>
#include "popkey.hpp"
#include <limits>
#include <cstdint>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include "config.h"

constexpr uint64_t HAPMAP_MNUM   = 3544454305642733928ULL;
constexpr uint64_t IHS2_JOB_MNUM = 521235431673ULL;
constexpr uint64_t IHS2_STD_MNUM = 5212354323174ULL;

#if defined(__GNUG__) && !defined(__ICC)
typedef unsigned long long v4ul __attribute__((vector_size(32)));
typedef unsigned long long v2ul __attribute__((vector_size(16)));
typedef unsigned int v8ui __attribute__((vector_size(32)));
typedef unsigned int v4ui __attribute__((vector_size(16)));
typedef unsigned short v16us __attribute__((vector_size(32)));
typedef unsigned short v8us __attribute__((vector_size(16)));
typedef short v16ss __attribute__((vector_size(32)));
typedef short v8ss __attribute__((vector_size(16)));
typedef int v8si __attribute__((vector_size(32)));
typedef char v32sc __attribute__((vector_size(32)));
typedef unsigned char v32uc __attribute__((vector_size(32)));
// typedef int v4si __attribute__((vector_size(16)));
#endif

//#define __AVX__
#ifdef VECLEN
#define VEC VECLEN
#else
#if defined(__AVX__)
#define VEC 4
#elif defined(__SSE2__)
#define VEC 2
#else
#define VEC 1
#endif
#endif

#if VEC==4
#define EQUAL(v1, v2) (v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2] && v1[3] == v2[3])
#define ZERO ((v4ul){0ULL, 0ULL, 0ULL, 0ULL})
#define BITSET_T_MAX ((v4ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__})
#define FIVES_VEC ((v4ul){0x5555555555555555,0x5555555555555555,0x5555555555555555,0x5555555555555555})
#define ONE_VEC ((v4ul){1,1,1,1})
#define VEC_ALIGN 32
#define POPCOUNT popcount4
#elif VEC==2
#define EQUAL(v1, v2) (v1[0] == v2[0] && v1[1] == v2[1])
#define ZERO ((v2ul){0ULL, 0ULL})
#define BITSET_T_MAX ((v2ul){__UINT64_MAX__,__UINT64_MAX__})
#define FIVES_VEC ((v2ul){0x5555555555555555,0x5555555555555555})
#define ONE_VEC ((v2ul){1,1})
#define VEC_ALIGN 16
#define POPCOUNT popcount2
#else
#define EQUAL(v1, v2) (v1 == v2)
#define ZERO 0ULL
#define BITSET_T_MAX __UINT64_MAX__
#define FIVES_VEC 0x5555555555555555
#define ONE_VEC 1
#define VEC_ALIGN 8
#define POPCOUNT popcount1
#endif

#if defined(__MINGW32__) && !defined(_ISOC11_SOURCE)
inline void* __attribute__((always_inline)) aligned_alloc(size_t alignment, size_t size)
{
    return __mingw_aligned_malloc(size, alignment);
}
inline void aligned_free(void* ptr)
{
    __mingw_aligned_free(ptr);
}
#elif (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !defined(_ISOC11_SOURCE)
inline void* __attribute__((always_inline)) aligned_alloc(size_t alignment, size_t size)
{
    void* ret;
    if (int r = posix_memalign(&ret, alignment, size) != 0)
    {
        throw std::bad_alloc();
    }
    return ret;
}
#define aligned_free free
#else
#define aligned_free free
#endif

template<typename T, std::size_t N>
constexpr std::size_t ctcBitsetSize()
{
    return (N/(sizeof(T)*8)) + ((N % (sizeof(T)*8) != 0) ? 1 : 0);
}

template<typename T>
std::size_t bitsetSize(std::size_t length)
{
    std::size_t bufferSize = length/(sizeof(T)*8);
    if (length % (sizeof(T)*8) != 0)
        bufferSize += 1;
    return bufferSize;
}

template<typename T, std::size_t N>
constexpr int ctcBitsetMaskShift()
{
    return (N == 0) ? 0 : (sizeof(T)*8 - (N % (sizeof(T)*8)));
}

template<typename T, std::size_t N>
constexpr T ctcBitsetMask()
{
    return (std::numeric_limits<T>::max() >> ctcBitsetMaskShift<T,N>());
}

template<typename T>
inline T bitsetMask(int length)
{
    return (std::numeric_limits<T>::max() >> (sizeof(T)*8-(length % (sizeof(T)*8))));
}

#ifdef __AVX__
inline v4ul bitsetMask4(int length)
{
    int len = length % 256;
    v4ul ret = (v4ul){0ULL,0ULL,0ULL,0ULL};
    if (len == 0)
    {
        ret = (v4ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__};
    }
    else if (len < 64)
    {
        ret[0] = bitsetMask<unsigned long long>(len);
    }
    else if (len < 128)
    {
        ret[0] = __UINT64_MAX__;
        ret[1] = bitsetMask<unsigned long long>(len);
    }
    else if (len < 192)
    {
        ret[0] = __UINT64_MAX__;
        ret[1] = __UINT64_MAX__;
        ret[2] = bitsetMask<unsigned long long>(len);
    }
    else
    {
        ret = (v4ul){__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__,__UINT64_MAX__};
        ret[3] = bitsetMask<unsigned long long>(len);
    }
    return ret;
}

#endif

#ifdef __SSE2__
inline v2ul bitsetMask2(int length)
{
    int len = length % 128;
    v2ul ret = (v2ul){0ULL,0ULL};
    if (len == 0)
    {
        ret = (v2ul){__UINT64_MAX__,__UINT64_MAX__};
    }
    else if (len < 64)
    {
        ret[0] = bitsetMask<unsigned long long>(len);
    }
    else
    {
        ret[0] = __UINT64_MAX__;
        ret[1] = bitsetMask<unsigned long long>(len);
    }
    return ret;
}
#endif

template<typename T, typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value>::type* = nullptr>
void convert(const char* line, T* buffer, std::size_t maxLength = 0)
{
    const std::size_t bits = sizeof(T)*8;

    std::size_t position = 0;
    for (std::size_t i = 0; line[i] != '\0'; ++i) {
        if (maxLength > 0 && position > maxLength)
            break;
        switch(line[i]) {
            case '0':
                position++;
                break;
            case '1':

                buffer[position/bits] |= ((T) (1ULL << (position % bits)));
                position++;
                break;
            case '\n':
            case '\r':
            case ' ':
                break;
            default:
                std::cout << "ERROR: Invalid character in ASCII haplotype map: " << line[i] << std::endl;
                throw std::runtime_error("ERROR: Not a valid ASCII haplotype map.");
                break;
        }
    }
}

template<typename T, typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value>::type* = nullptr>
void convert(const char* line, T* buffer, const PopKey& pk, std::size_t maxLength = 0)
{
    const std::size_t bits = sizeof(T)*8;
    std::size_t position = 0;
    std::size_t pk_pos = 0;
    for (std::size_t i = 0; line[i] != '\0'; ++i) {
        if (maxLength > 0 && position > maxLength)
            break;
        if (position % bits == 0)
            assert(buffer[position/bits] == 0ULL);
        switch(line[i]) {
            case '0':
                if (pk[i/2])
                {
                    buffer[position/bits] &= ((T) -1) ^ ((T) (1ULL << (position % bits)));
                    position++;
                }
                ++pk_pos;
                break;
            case '1':
                if (pk[i/2])
                {
                    buffer[position/bits] |= ((T) (1ULL << (position % bits)));
                    position++;
                }
                ++pk_pos;
                break;
            case '\n':
            case '\r':
            case ' ':
                break;
            default:
                std::cout << "ERROR: Invalid character in ASCII haplotype map: " << line[i] << std::endl;
                throw std::runtime_error("ERROR: Not a valid ASCII haplotype map.");
                break;
        }
    }
}

struct Stats
{
    Stats() : mean{}, stddev{} {}
    double mean;
    double stddev;
};

template<typename T, typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value>::type* = nullptr>
void print_snp(T* buffer, std::size_t length)
{
    constexpr std::size_t bits = sizeof(T)*8;
    for (std::size_t i = 0; i < length; ++i)
    {
        if (buffer[i/bits] & (1ULL << (i%bits)))
            std::cout << "1 ";
        else
            std::cout << "0 ";
    }
    std::cout << std::endl;
}

void filter(unsigned long long* in, unsigned long long* out, std::size_t inlen, std::size_t outlen, const PopKey& pk);
double binom_2(double n);
double nearest(double target, double number);
Stats stats(const std::vector<double>& list);
std::vector<std::string> splitString(const std::string input, char delim);

inline int popcount1(unsigned long long val)
{
    return __builtin_popcountll(val);
}

inline int popcount4(v4ul val)
{
    int ret = __builtin_popcountll(val[0]);
    ret += __builtin_popcountll(val[1]);
    ret += __builtin_popcountll(val[2]);
    ret += __builtin_popcountll(val[3]);
    return ret;
}

inline int popcount2(v2ul val)
{
    int ret = __builtin_popcountll(val[0]);
    ret += __builtin_popcountll(val[1]);
    return ret;
}

#endif // HAPBIN_HPP

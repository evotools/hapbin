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

#ifndef COMMON_HPP
#define COMMON_HPP

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__  << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(ERR_ASSERT); \
        } \
    } while (false)
#else /* NDEBUG */
#define ASSERT(condition, message) do {} while(false)
#endif

namespace mpirpc
{

/**
 * Passer can be used along with uniform initialization to unpack parameter packs
 * and execute the parameters in the order in which they appear. This is necessary
 * for correctness when side effects are important.
 */
struct Passer {
    Passer(...) {}
};

template<typename F>
struct FunctionParts;

template<typename R, class Class, typename... Args>
struct FunctionParts<R(Class::*)(Args...)>
{
    using return_type = R;
    using class_type = Class;
    using function_type = R(Class::*)(Args...);
};

template<typename R, typename... Args>
struct FunctionParts<R(*)(Args...)>
{
    using return_type = R;
    using function_type = R(*)(Args...);
};

using FunctionHandle = unsigned long long;
using TypeId = unsigned long long;
using ObjectId = unsigned long long;

}

#endif // COMMON_HPP

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

#ifndef ORDEREDCALL_HPP
#define ORDEREDCALL_HPP

#include "common.hpp"

namespace mpirpc {

/**
 * @brief The OrderedCall<F> class
 *
 * OrderedCall<F> can be used to call functions such that the
 * side effects of the parameters are evaluated in the order
 * in which they appear. This is done using uniform initilization:
 * constructing using {} brackets instead of (). This bound call
 * can then be executed using OrderedCall<F>::operator().
 */
template <typename F>
struct OrderedCall;

/**
 * Specialization of OrderedCall<F> for function pointer calls.
 */
template<typename R, typename... Args>
struct OrderedCall<R(*)(Args...)>
{
    OrderedCall(R(*function)(Args...), Args&&... args)
    {
        bound = std::bind(function, std::forward<Args>(args)...);
    }

    R operator()()
    {
        return bound();
    }

    std::function<R()> bound;
};

/**
 * Specialization of OrderedCall<F> for member function pointer calls.
 */
template<typename R, typename Class, typename... Args>
struct OrderedCall<R(Class::*)(Args...)>
{
    OrderedCall(R(Class::*function)(Args...), Class *c, Args&&... args)
    {
        bound = std::bind(function, c, std::forward<Args>(args)...);
    }

    R operator()()
    {
        return bound();
    }

    std::function<R()> bound;
};

/**
 * Specialization of OrderedCall<F> for std::function calls.
 */
template<typename R, typename... Args>
struct OrderedCall<std::function<R(Args...)>>
{
    OrderedCall(std::function<R(Args...)> &function, Args&&... args)
    {
        bound = std::bind(function, std::forward<Args>(args)...);
    }

    R operator()()
    {
        return bound();
    }

    std::function<R()> bound;

};

}

#endif // ORDEREDCALL_HPP

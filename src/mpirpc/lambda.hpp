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

#ifndef LAMBDA_HPP
#define LAMBDA_HPP

#include <functional>

namespace mpirpc {

template <typename F>
struct LambdaTraits : public LambdaTraits<decltype(&F::operator())>
{};

template <typename C, typename R, typename... Args>
struct LambdaTraits<R(C::*)(Args...) const>
{
    //using lambda_fnPtr       = R(*)(Args...);
    using lambda_stdfunction = std::function<R(Args...)>;
};

}

#endif // LAMBDA_HPP

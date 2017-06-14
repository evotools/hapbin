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
 
#ifndef MPITYPE_H
#define MPITYPE_H

#include <mpi.h>
#include <cstdint>

template<typename T>
MPI_Datatype mpiType();

template<>
MPI_Datatype mpiType<double>();

template<>
MPI_Datatype mpiType<long double>();

template<>
MPI_Datatype mpiType<float>();

template<>
MPI_Datatype mpiType<char>();

template<>
MPI_Datatype mpiType<unsigned char>();

template<>
MPI_Datatype mpiType<short>();

template<>
MPI_Datatype mpiType<unsigned short>();

template<>
MPI_Datatype mpiType<int>();

template<>
MPI_Datatype mpiType<unsigned int>();

template<>
MPI_Datatype mpiType<long>();

template<>
MPI_Datatype mpiType<unsigned long>();

template<>
MPI_Datatype mpiType<long long>();

template<>
MPI_Datatype mpiType<unsigned long long>();

#endif /* MPITYPE_H */

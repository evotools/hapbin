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

#ifndef REDUCE_HPP
#define REDUCE_HPP

#include <type_traits>
#include <mpi.h>

namespace mpirpc
{

char reduce(char value, MPI_Op op, int count, int root, MPI_Comm comm);
short reduce(short value, MPI_Op op, int count, int root, MPI_Comm comm);
unsigned short reduce(unsigned short value, MPI_Op op, int count, int root, MPI_Comm comm);
int reduce(int value, MPI_Op op, int count, int root, MPI_Comm comm);
unsigned int reduce(unsigned int value, MPI_Op op, int count, int root, MPI_Comm comm);
long reduce(long value, MPI_Op op, int count, int root, MPI_Comm comm);
unsigned long reduce(unsigned long value, MPI_Op op, int count, int root, MPI_Comm comm);
float reduce(float value, MPI_Op op, int count, int root, MPI_Comm comm);
double reduce(double value, MPI_Op op, int count, int root, MPI_Comm comm);

char allreduce(char value, MPI_Op op, int count, MPI_Comm comm);
short allreduce(short value, MPI_Op op, int count, MPI_Comm comm);
unsigned short allreduce(unsigned short value, MPI_Op op, int count, MPI_Comm comm);
int allreduce(int value, MPI_Op op, int count, MPI_Comm comm);
unsigned int allreduce(unsigned int value, MPI_Op op, int count, MPI_Comm comm);
long allreduce(long value, MPI_Op op, int count, MPI_Comm comm);
unsigned long allreduce(unsigned long value, MPI_Op op, int count, MPI_Comm comm);
float allreduce(float value, MPI_Op op, int count, MPI_Comm comm);
double allreduce(double value, MPI_Op op, int count, MPI_Comm comm);

}

#endif // REDUCE_HPP

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

#include "reduce.hpp"

namespace mpirpc {

char reduce(char value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    char ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_BYTE, op,  root, comm);
    return ret;
}

short reduce(short value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    short ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_SHORT, op,  root, comm);
    return ret;
}

unsigned short reduce(unsigned short value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    unsigned short ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_UNSIGNED_SHORT, op,  root, comm);
    return ret;
}

int reduce(int value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    int ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_INT, op,  root, comm);
    return ret;
}

unsigned int reduce(unsigned int value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    unsigned int ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_UNSIGNED, op,  root, comm);
    return ret;
}

long reduce(long value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    long ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_LONG, op,  root, comm);
    return ret;
}

unsigned long reduce(unsigned long value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    unsigned long ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_UNSIGNED_LONG, op,  root, comm);
    return ret;
}

float reduce(float value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    float ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_FLOAT, op,  root, comm);
    return ret;
}

double reduce(double value, MPI_Op op, int count, int root, MPI_Comm comm)
{
    double ret = 0;
    MPI_Reduce((void*) &value, &ret, count, MPI_DOUBLE, op,  root, comm);
    return ret;
}

///////////////////////

char allreduce(char value, MPI_Op op, int count, MPI_Comm comm)
{
    char ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_BYTE, op, comm);
    return ret;
}

short allreduce(short value, MPI_Op op, int count, MPI_Comm comm)
{
    short ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_SHORT, op, comm);
    return ret;
}

unsigned short allreduce(unsigned short value, MPI_Op op, int count, MPI_Comm comm)
{
    unsigned short ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_UNSIGNED_SHORT, op, comm);
    return ret;
}

int allreduce(int value, MPI_Op op, int count, MPI_Comm comm)
{
    int ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_INT, op, comm);
    return ret;
}

unsigned int allreduce(unsigned int value, MPI_Op op, int count, MPI_Comm comm)
{
    unsigned int ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_UNSIGNED, op, comm);
    return ret;
}

long allreduce(long value, MPI_Op op, int count, MPI_Comm comm)
{
    long ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_LONG, op, comm);
    return ret;
}

unsigned long allreduce(unsigned long value, MPI_Op op, int count, MPI_Comm comm)
{
    unsigned long ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_UNSIGNED_LONG, op, comm);
    return ret;
}

float allreduce(float value, MPI_Op op, int count, MPI_Comm comm)
{
    float ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_FLOAT, op, comm);
    return ret;
}

double allreduce(double value, MPI_Op op, int count, MPI_Comm comm)
{
    double ret = 0;
    MPI_Allreduce((void*) &value, &ret, count, MPI_DOUBLE, op, comm);
    return ret;
}



}

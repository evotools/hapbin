# MPIRPC: MPI based invocation of functions on other ranks
# Copyright (C) 2014-2017 Colin MacLean <cmaclean@illinois.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
set(SRC_LIST manager.cpp objectwrapper.cpp parameterstream.cpp mpitype.cpp)
add_library(mpirpc STATIC ${SRC_LIST})
target_link_libraries(mpirpc ${MPI_CXX_LIBRARIES})

install(TARGETS mpirpc DESTINATION lib)
install(FILES common.hpp lambda.hpp manager.hpp objectwrapper.hpp orderedcall.hpp parameterstream.hpp mpitype.hpp DESTINATION include/mpirpc)

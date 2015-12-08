# Hapbin: A fast binary implementation EHH, iHS, and XPEHH
# Copyright (C) 2014  Colin MacLean <s0838159@sms.ed.ac.uk>
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

project(hapbin)
cmake_minimum_required(VERSION 2.6)
option(MARCH_NATIVE "Optimize with -march=native" ON)
option(USE_MPI "Build with MPI, if available" ON)
option(USE_LTO "Build with Link Time Optimizations" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(GetGitRevisionDescription)
git_describe(VERSION --tags --dirty=-dirty)

#parse the version information into pieces.
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH
    "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+(.*)" "\\1" VERSION_SHA1
    "${VERSION}")
set(VERSION_SHORT "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
if(CMAKE_COMPILER_IS_GNUCXX)
    if(USE_LTO)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto -fuse-linker-plugin -fuse-ld=gold")
        set(CMAKE_AR "gcc-ar")
        set(CMAKE_NM "gcc-nm")
        set(CMAKE_RANLIB "gcc-ranlib")
    endif(USE_LTO)
    if(MARCH_NATIVE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif(MARCH_NATIVE)
endif(CMAKE_COMPILER_IS_GNUCXX)

find_package(OpenMP)

#add_definitions(-std=c++11 -Wall -O3 -march=corei7)

include_directories("${PROJECT_BINARY_DIR}")

include(FindMPI)
if(MPI_FOUND AND USE_MPI)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    include_directories("${MPI_CXX_INCLUDE_PATH}")
else(MPI_FOUND AND USE_MPI)
    set(MPI_C_LIBRARIES "")
    set(MPI_CXX_LIBRARIES "")
endif(MPI_FOUND AND USE_MPI)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${OpenMP_CXX_FLAGS}")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")

if(MPI_FOUND AND USE_MPI)
   set(HAVE_MPILIB 1)
   set(MPIRPC_FOUND 1)
   add_subdirectory("mpirpc")
   include_directories("mpirpc")
endif(MPI_FOUND AND USE_MPI)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")

set(core_SRCS ehhfinder.cpp ihsfinder.cpp hapmap.cpp hapbin.cpp)
add_library(hapbin SHARED ${core_SRCS})
set_target_properties(hapbin PROPERTIES VERSION 0 SOVERSION 0.0.0)

set(ihsbin_SRCS main-ihs.cpp)
add_executable(ihsbin ${ihsbin_SRCS})

set(ehhbin_SRCS main-ehh.cpp)
add_executable(ehhbin ${ehhbin_SRCS})

set(xpehhbin_SRCS main-xpehh.cpp)
add_executable(xpehhbin ${xpehhbin_SRCS})

set(hapbinconv_SRCS main-conv.cpp)
add_executable(hapbinconv ${hapbinconv_SRCS})

if(MPI_FOUND AND USE_MPI)
    target_link_libraries(ihsbin hapbin mpirpc ${MPI_CXX_LIBRARIES})
    target_link_libraries(xpehhbin hapbin mpirpc ${MPI_CXX_LIBRARIES})
else(MPI_FOUND AND USE_MPI)
    target_link_libraries(ihsbin hapbin)
    target_link_libraries(xpehhbin hapbin)
endif(MPI_FOUND AND USE_MPI)
target_link_libraries(ehhbin hapbin)
target_link_libraries(hapbinconv hapbin)

install(TARGETS hapbin DESTINATION lib)
install(TARGETS ihsbin DESTINATION bin)
install(TARGETS ehhbin DESTINATION bin)
install(TARGETS xpehhbin DESTINATION bin)
install(TARGETS hapbinconv DESTINATION bin)
install(FILES calcmpiselect.hpp calcnompiselect.hpp calcselect.hpp argparse.hpp hapmap.hpp hapbin.hpp ihsfinder.hpp ehhfinder.hpp DESTINATION include/hapbin)

include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Hapbin is a fast and efficient implementation of EHH and iHS calculations using a bitwise algorithm.")
set(CPACK_PACKAGE_VENDOR "Colin MacLean, University of Edinburgh")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "x86_64")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Colin MacLean <s0838159@sms.ed.ac.uk>")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "mpich (>= 3.0.0), libgcc1 (>= 1:4.7), libstdc++6 (>= 4.7)")
set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
set(CPACK_RPM_PACKAGE_REQUIRES "mpich >= 3.0.0, libgcc1 >= 4.7, libstdc++6 >= 4.7") 
include(CPack)

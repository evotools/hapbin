# Hapbin: A fast binary implementation EHH, iHS, and XPEHH
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

project(hapbin)
cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(MARCH_NATIVE "Optimize with -march=native" ON)
option(MAKE_TESTS "Build tests" OFF)
option(USE_MPI "Build with MPI, if available" ON)
option(USE_LTO "Build with Link Time Optimizations" OFF)
#set(CMAKE_CXX_COMPILER "clang++")
#set(CMAKE_CXX_COMPILER "g++")
set(ENV{CRAYPE_LINK_TYPE} "dyanamic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")

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

macro(USE_CXX11)
  if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
  else ()
    set (CMAKE_CXX_STANDARD 11)
  endif ()
endmacro(USE_CXX11)

USE_CXX11()

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")

if(MPI_FOUND AND USE_MPI)
   set(HAVE_MPILIB 1)
   set(MPIRPC_FOUND 1)
   add_subdirectory("mpirpc")
   include_directories("mpirpc")
endif(MPI_FOUND AND USE_MPI)

configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")
set(core_SRCS ehhfinder.cpp ihsfinder.cpp hapmap.cpp hapbin.cpp ehhpairfinder.cpp ihspairfinder.cpp ihspairjob.cpp popkey.cpp calcselect.hpp calcmpiselect.hpp)
add_library(hapbin SHARED ${core_SRCS})
set_target_properties(hapbin PROPERTIES VERSION 0 SOVERSION 0.0.0)

set(ihsbin_SRCS main-ihs.cpp)
add_executable(ihsbin ${ihsbin_SRCS})

set(ehhbin_SRCS main-ehh.cpp)
add_executable(ehhbin ${ehhbin_SRCS})

set(ehh2bin_SRCS main-ehh2.cpp)
add_executable(ehh2bin ${ehh2bin_SRCS})

set(xpehhbin_SRCS main-xpehh.cpp)
add_executable(xpehhbin ${xpehhbin_SRCS})

set(hapbinconv_SRCS main-conv.cpp)
add_executable(hapbinconv ${hapbinconv_SRCS})

set(ihs2binsub_SRCS main-ihs2sub.cpp)
add_executable(ihs2binsub ${ihs2binsub_SRCS})

if(MPI_FOUND AND USE_MPI)
    set(ihs2bin_SRCS main-ihs2.cpp)
    add_executable(ihs2bin ${ihs2bin_SRCS})
endif(MPI_FOUND AND USE_MPI)

if(MPI_FOUND AND USE_MPI)
    target_link_libraries(ihsbin hapbin mpirpc ${MPI_CXX_LIBRARIES})
    target_link_libraries(ihs2bin hapbin mpirpc ${MPI_CXX_LIBRARIES})
    target_link_libraries(xpehhbin hapbin mpirpc ${MPI_CXX_LIBRARIES})
else(MPI_FOUND AND USE_MPI)
    target_link_libraries(ihsbin hapbin)
    target_link_libraries(xpehhbin hapbin)
endif(MPI_FOUND AND USE_MPI)
target_link_libraries(ehhbin hapbin)
target_link_libraries(ehh2bin hapbin)
target_link_libraries(hapbinconv hapbin)
target_link_libraries(ihs2binsub hapbin)

if (MAKE_TESTS)
#    add_executable(hapmap_test hapmap_test.cpp)
#    target_link_libraries(hapmap_test hapbin CppUTest)
#    add_test(NAME hapmap_test COMMAND hapmap_test)
endif (MAKE_TESTS)

install(TARGETS hapbin DESTINATION lib)
install(TARGETS ihsbin DESTINATION bin)
install(TARGETS ehhbin DESTINATION bin)
install(TARGETS ehh2bin DESTINATION bin)
install(TARGETS xpehhbin DESTINATION bin)
install(TARGETS hapbinconv DESTINATION bin)
install(TARGETS ihs2binsub DESTINATION bin)
if(USE_MPI AND MPI_FOUND)
    install(TARGETS ihs2bin DESTINATION bin)
endif(USE_MPI AND MPI_FOUND)
install(FILES calcmpiselect.hpp calcnompiselect.hpp calcselect.hpp
    argparse.hpp hapmap.hpp hapbin.hpp ihsfinder.hpp ehhfinder.hpp
    ihspairfinder.hpp ehhpairfinder.hpp DESTINATION include/hapbin)

include(InstallRequiredSystemLibraries)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Hapbin is a fast and efficient implementation of EHH and iHS calculations using a bitwise algorithm.")
set(CPACK_PACKAGE_VENDOR "Colin MacLean, University of Edinburgh")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/../README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "2")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "x86_64")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Colin MacLean <cmaclean@illinois.edu>")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "mpich (>= 3.0.0), libgcc1 (>= 1:4.7), libstdc++6 (>= 4.7)")
set(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
set(CPACK_RPM_PACKAGE_REQUIRES "mpich >= 3.0.0, libgcc1 >= 4.7, libstdc++6 >= 4.7")
include(CPack)

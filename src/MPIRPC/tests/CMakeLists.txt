include_directories(${Qt5Core_INCLUDE_DIRS} ${Qt5Test_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR})
set(CMAKE_EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

set(mpirpctest_SRCS mpirpctest.cpp ../manager.cpp ../manager.hpp ../common.hpp ../lambda.hpp
    ../objectwrapper.hpp ../objectwrapper.cpp ../orderedcall.hpp ../reduce.hpp ../reduce.cpp
    ../parameterstream.cpp ../parameterstream.hpp)
add_executable(mpirpcTest ${mpirpctest_SRCS})
add_test(mpirpcTest mpirpcTest)
target_link_libraries(mpirpcTest ${Qt5Core_LIBRARIES} ${Qt5Test_LIBRARIES} ${MPI_CXX_LIBRARIES})
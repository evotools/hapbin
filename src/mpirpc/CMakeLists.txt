set(SRC_LIST manager.cpp objectwrapper.cpp parameterstream.cpp mpitype.cpp)
add_library(mpirpc STATIC ${SRC_LIST})
target_link_libraries(mpirpc ${MPI_CXX_LIBRARIES})

install(TARGETS mpirpc DESTINATION lib)
install(FILES common.hpp lambda.hpp manager.hpp objectwrapper.hpp orderedcall.hpp parameterstream.hpp mpitype.hpp DESTINATION include/mpirpc)

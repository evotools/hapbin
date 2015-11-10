module switch PrgEnv-cray PrgEnv-gnu
module load gcc
module load cmake
cmake ../src/ -DCMAKE_TOOLCHAIN_FILE=BW-toolchain.cmake #-DUSE_MPI=OFF
make -j10

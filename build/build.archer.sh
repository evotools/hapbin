module load gcc
module load cmake
module switch PrgEnv-cray PrgEnv-gnu
cmake ../src/ -DCMAKE_TOOLCHAIN_FILE=ARCHER-toolchain.cmake
make -j10

hapbin
======

`hapbin` is a collection of tools for efficiently calculating [Extended Haplotype Homozygosity (EHH)](http://dx.doi.org/10.1038/nature01140), the [Integrated Haplotype Score (iHS)](http://dx.doi.org/10.1371/journal.pbio.0040072) and the [Cross Population Extended Haplotype Homozogysity (XP-EHH)](http://www.nature.com/nature/journal/v449/n7164/full/nature06250.html) statistic.

## Tools ##

The `hapbin` suite contains the following tools:

   * `ehhbin --hap [.hap/.hapbin file] --map [.map file] --locus [locus] --out [output prefix]` - calculate the EHH
   * `ihsbin --hap [.hap/.hapbin file] --map [.map file] --out [output prefix]` - calculate the iHS of all loci in a `.hap/.hapbin` file
   * `xpehhbin --hapA [Population A .hap/.hapbin] --hapB [Population B .hap/.hapbin] --map [.map file] --out [output prefix]` - calculate the XPEHH of all loci in `.hap/.hapbin` files.
   * `hapbinconv --hap [.hap ASCII file] --out [.hapbin binary file]` - convert .hap file to more size efficient binary format.

For additional options, see `[executable] --help`.

## Copyright and License ##

This code is licensed under the GPL v3. Copyright is retained by the original authors, Colin Maclean and the University of Edinburgh.

## Building from source code ##

### Dependencies ###

   * A C++11 capable compiler (GCC >= 4.7 for required features). OpenMP support required for threaded execution.

   * Optional dependency: MPI for execution on distributed memory systems (clusters/supercomputers).

If any of these are not already installed on your system then for the main Linux distributions they can simply be added via their package managers.

For example to install on **Ubuntu** (tested on 14.04 LTS):

     sudo apt-get update
     sudo apt-get install git cmake libcr-dev mpich libmpich-dev

On **openSUSE** (tested on Enterprise Server 12):

     sudo zypper install cmake git-core gcc-c++ openmpi openmpi-devel
     export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/lib64:/usr/lib64/mpi/gcc/openmpi/bin

On **Red Hat** (tested on Enterprise Linux 7.1):

     sudo yum install cmake git gcc-c++

For those servers where CMake is not installed, and you do not have the necessary permissions to add it as above, precompiled binaries can be downloaded from [here](http://www.cmake.org/download/).

### Building the source code ###

An out of source build is suggested in order to keep the source directory clean. To do this, check out the hapbin source, move to the build directory, then run `cmake [path to src directory]`. Once CMake has finished generating the necessary files, simply run `make`.

For example:

     git clone https://github.com/evotools/hapbin.git
     cd hapbin/build/
     cmake ../src/
     make

### Installing with an alternate toolchain ###

First, check out the hapbin source:
    
    git clone https://github.com/evotools/hapbin.git
    cd hapbin/build/

Next, create a `toolchain.cmake` file with the necessary overrides:

    ...
    SET(CMAKE_C_COMPILER "/path/to/c/compiler")
    SET(CMAKE_CXX_COMPILER "/path/to/cxx/compiler")
    ...

`MPI_C_LIBRARIES`, `MPI_CXX_LIBRARIES`, `MPI_C_INCLUDE_PATH`, and `MPI_CXX_INCLUDE_PATH` can be set in this file, too, if necessary.

Then, tell cmake to use this toolchain and build:

    cmake ../src/ -DCMAKE_TOOLCHAIN_FILE=toolchain.cmake
    make

### Installing on ARCHER ###

If you are using hapbin on the [ARCHER UK National HPC Service](http://www.archer.ac.uk/), follow these steps:

   1. Download hapbin: `git clone https://github.com/evotools/hapbin.git`

   2. Navigate to the build directory: `cd hapbin/build`

   3. Run `. build.archer.sh` or `source build.archer.sh` to load/switch required environment modules and configure/install `hapbin`.

   4. Install to `/work/[project code]/[group code]/[username]/hapbin/`  by typing `make install`

   5. Copy desired `hapbin.[haps].[threads].pbs` from `hapbin/tools/pbs/hapbin` to `/work/[project code]/[group code]/[username]/hapbin/bin/`

   6. `cd /work/[project code]/[group code]/[username]/hapbin/bin`, edit PBS as desired, and submit to the batch queue with `qsub`.

### Input file formats ###

The hap files (`--hap`), containing phased haplotypes, should be in IMPUTE [hap format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-h). These can be optionally converted to smaller binary files for use with the hapbin suite of tools using `hapbinconv`. IMPUTE provides phased haplotypes in this format for several publically available human cohorts [here](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference). If your data is in VCF format it can be converted to IMPUTE format using [vcftools](https://vcftools.github.io).

The map files (`--map`) should be in the same format as used by [Selscan](https://github.com/szpiech/selscan) with one row per variant and four space-separated columns specifiying chromosome, locus ID, genetic position and physical position.

### Output file formats ###

- ehhbin outputs five columns. The first three being the locus' ID and its genetic and physical positions. These are followed by two columns corresponding to the EHH for each of the alleles at this locus (allele coded as 0 then 1).
- ihsbin outputs two files, the first containing unstandardised iHS for allele 0 and the second (with the .std extension) containing the corresponding standardised iHS (alleles grouped in to 2% frequency bins for standardisation by default). Each of these output files contains two columns: the SNP locus id (as specified in the map file) and corresponding iHS value.
- xpehh also outputs a file containing two columns: the SNP locus id (as specified in the map file) and corresponding XP-EHH value.

### Examples ###

Example command for calculating EHH for a variant with ID (`--locus`) of 140465 as specified in the map input file. Output is redirected to a file named 140465_EHH.txt:

     ./ehhbin --hap 1000GP_Phase3.GBR.chr22.hap --map chr22.map --locus 140465 > 140465_EHH.txt

Example command for calculating the iHS of all variants with a minor allele frequency greater than 10% (`--minmaf 0.1`) and specifying that the integral of the observed decay of EHH (i.e. iHH, see [Voight et al.](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072) for more information) should be calculated up to the point at which EHH drops below 0.1 (`--cutoff 0.1`):

     ./ihsbin --hap 1000GP_Phase3.GBR.chr22.hap --map chr22.map --out chr22_iHS --minmaf 0.1 --cutoff 0.1

Example command for calculating XP-EHH with default values for minor allele frequency and EHH cutoff:

     ./xpehhbin --hapA 1000GP_Phase3.GBR.chr22.hap --hapB 1000GP_Phase3.YRI.chr22.hap --map chr22.map --out chr22_GBRvsYRI_XPEHH

Each of the input files referred to in these examples can be found in the data directory.

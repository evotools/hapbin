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

### Input file formats ###

The hap files (`--hap`), containing phased haplotypes, should be in IMPUTE [hap format](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#-h). These can be optionally converted to smaller binary files for use with the hapbin suite of tools using hapbincov. IMPUTE provides phased haplotypes in this format for several publically available human cohorts [here](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference).

The map files (`--map`) should be in the same format as used by [Selscan](https://github.com/szpiech/selscan) with one row per variant and four space-separated columns specifiying chromosome, locus ID, genetic position and physical position.

### Output file formats ###

- ehhbin outputs two columns, the EHH for each allele (0 and 1) at each location.
- ihsbin outputs two files, the first containing unstandardised iHS for allele 0 and the second (with the .std extension) containing the corresponding standardised iHS (alleles grouped in to 2% frequency bins for standardisation by default). Each of these output files contains two columns: the SNP locus id (as specified in the map file) and corresponding iHS value.
- xpehh output file also contains two columns: the SNP locus id (as specified in the map file) and corresponding XP-EHH value.

### Examples ###

Example command for calculating EHH for a variant with ID (`--locus`) of 9189 as specified in the map input file. Output is redirected to file named 9189_EHH.txt:

     ehhbin --hap phasedHaplotypes_chr22.hap --map chr22.map --locus 9189 > 9189_EHH.txt

Example command for calculating the iHS of all variants with a minor allele frequency greater than 10% (`--minmaf 0.1`) and specifying that the integral of the observed decay of EHH (i.e. iHH, see [Voight et al.](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0040072) for more information) should be calculated up to the point at which EHH drops below 0.1 (`--cutoff 0.1`):

     ihsbin --hap phasedHaplotypes_chr22.hap --map chr22.map --out chr22_iHS --minmaf 0.1 --cutoff 0.1

Example command for calculating XP-EHH with default values for minor allele frequency and EHH cutoff:

     xpehhbin --hapA EUR_phasedHaplotypes_chr22.hap --hapB AFR_phasedHaplotypes_chr22.hap --map chr22.map --out chr22_EURvsAFR_XPEHH


## Copyright and License ##

This code is licensed under the GPL v3. Copyright is retained by the original authors, Colin Maclean and the University of Edinburgh.

## Building from source code ##

### Dependencies ###

   * A C++11 capable compiler (GCC >= 4.7 for required features). OpenMP support required for threaded execution. Requires MPIRPC from https://github.com/camaclean/MPIRPC.

   * Optional dependency: MPI for execution on distributed memory systems (clusters/supercomputers).

   * Optional dependency: QT 5 for the unit test framework. On Ubuntu, see: https://qt-project.org/wiki/Install_Qt_5_on_Ubuntu.

### Building the source code ###

An out of source build is suggested in order to keep the source directory clean. To do this, create a build directory, then run `cmake [path to directory]` in the build directory.

For example:

     cd /path/to/hapbin
     cd build
     cmake ../src/

Once CMake has finished generating the necessary files, simply run `make`.

The test programs are created in a `test` subdirectory. Run these test programs with `-help` or see the Qt 5 QTest framework documentation for testing and benchmarking options.

Running `ctest` or `make test` will run all test programs.

### Installing on Ubuntu ###

First ensure packages required for obtaining and compiling code are installed as well as mpi packages used for parallelisation if required.

     sudo apt-get update
     sudo apt-get install git cmake libcr-dev mpich libmpich-dev

Install [MPIRPC](https://github.com/camaclean/MPIRPC) to chosen directory.

     git clone https://github.com/camaclean/MPIRPC.git
     cd MPIRPC/build/
     cmake ../src/
     make
     sudo make install

Finally download and compile hapbin.

     cd ../../
     git clone -b master https://github.com/evotools/hapbin.git
     cd hapbin/build/
     cmake ../src/
     make -j 4

#### Installing without root permissions ####

If you dont want or unable to install MPIRPC system-wide, for example if you do not have the required permissions, you can specify the directory for it to be installed into and then direct cmake to this directory when compiling hapbin. For example to install MPIRPC to the subdirectory MPIRPC/install in your home directory:

     git clone https://github.com/camaclean/MPIRPC.git
     cd MPIRPC/build/
     cmake -DCMAKE_INSTALL_PREFIX=$HOME/MPIRPC/install/ ../src/
     make
     make install
     
     cd ../../
     git clone -b master https://github.com/evotools/hapbin.git
     cd hapbin/build/
     cmake -DCMAKE_PREFIX_PATH=$HOME/MPIRPC/install/ ../src/
     make -j 4

The two instances of $HOME/MPIRPC/install/ in the above code can be changed to whichever directory you prefer.

### Installing on ARCHER ###

If you are using hapbin on the [ARCHER UK National HPC Service](http://www.archer.ac.uk/), follow these steps:

   1. Install [MPIRPC](https://github.com/camaclean/MPIRPC).

   2. Navigate to the build directory: `cd hapbin/build`

   3. Run `. build.archer.sh` or `source build.archer.sh` to load/switch required environment modules and configure/install `hapbin`.

   4. Install to `/work/[project code]/[group code]/[username]/hapbin/`  by typing `make install`

   5. Copy desired `hapbin.[haps].[threads].pbs` from `hapbin/tools/pbs/hapbin` to `/work/[project code]/[group code]/[username]/hapbin/bin/`

   6. `cd /work/[project code]/[group code]/[username]/hapbin/bin`, edit PBS as desired, and submit to the batch queue with `qsub`.

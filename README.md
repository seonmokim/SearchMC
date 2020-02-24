# SearchMC

## Overview

SearchMC is an approximate model counter using XOR streamlining techniques. SearchMC takes an SMT formula or a CNF formula as input and gives a lower and a upper bound of the true model influence (log base 2 of the true model counts) with a given confidence level. 

## Paper

* [Bit-Vector Model Counting using Statistical Estimation (TACAS 2018)](https://link.springer.com/chapter/10.1007/978-3-319-89960-2_8/)

## Install SearchMC

* Build particle-filter.c
```
$ gcc asa241.c particle-filter.c -o particle-filter -lm
```
This binary file (``particle-filter``) should be in SearchMC folder.

* Install STP

SearchMC uses STP for converting an SMT formula to a CNF formula. STP relies on a few dependencies, namely boost, flex/bison and minisat. Installing the required header files and binaries can be achieved through the following. (Tested on Ubuntu 14.04.)
```
$ sudo apt-get install cmake bison flex libboost-all-dev python perl
```
Installing minisat can be achieved by running
```
$ git clone https://github.com/niklasso/minisat.git
$ cd minisat
$ mkdir build && cd build
$ cmake ..
$ make
$ sudo make install
$ sudo ldconfig
```
Installing STP with printing out corresponding output variables from an SMT formula can be achieved by running
```
$ git clone https://github.com/seonmokim/stp.git
$ cd stp
$ mkdir build && cd build
$ cmake -G 'Unix Makefiles' /path/to/stp/source/root
$ make
```
Then, copy the STP binary file (``stp-2.1.2``) into SearchMC folder.

* Install Cryptominisat2 with projection

Make sure to build with Gaussian Eliminiation.
```
$ git clone https://github.com/seonmokim/cryptominisat2.git
$ cd cryptominisat
$ mkdir build && cd build
$ ../configure
$ make
```
Then, copy the Cryptominisat binary file (``cryptominisat``) into SearchMC folder.

* Install Z3

Make sure to build with Gaussian Eliminiation.
```
$ git clone https://github.com/Z3Prover/z3.git
$ cd z3
$ python scripts/mk_make.py
$ cd build
$ make
$ sudo make install
```
Then, copy the Z3 binary file (``z3``) into SearchMC folder.

* Install MathSAT5

MathSAT5 is available at http://mathsat.fbk.eu/

Copy the MathSAT5 binary file (``mathsat``) into SearchMC folder.

* In-place version of build instructions

It is also possible to compile all of the prerequisite software
in-place if you cannot or do not wish to install anything into your
system directories. This requires extra options to tell the programs
about the locations of their prerequisites, and it works most easily
if everything is compiled using static linking. Below are some sample
command lines, for building everything in a subdirectory ``builds`` of
the directory where these instructions are found:

```
$ mkdir builds
$ cd builds
$ git clone https://github.com/niklasso/minisat.git
$ cd minisat
$ mkdir build install
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=../install
$ make
$ make install
$ cd ../..
$ git clone https://github.com/seonmokim/stp.git
$ cd stp
$ mkdir build
$ cd build
$ cmake -DMINISAT_LIBRARY=../../minisat/install/lib/libminisat.a -DMINISAT_INCLUDE_DIR=../../minisat/install/include -DBUILD_STATIC_BIN=ON -DBUILD_SHARED_LIBS=OFF -DENABLE_PYTHON_INTERFACE=OFF ..
$ make
$ cp stp-2.1.2 ../../..
$ cd ../..
$ git clone https://bitbucket.org/malb/m4ri.git
$ cd m4ri
$ mkdir install
$ autoreconf --install
$ ./configure --prefix=$(pwd)/install
$ make
$ make install
$ cd ..
$ git clone https://gitlab.com/QIF/cryptominisat4.git cryptominisat5
$ cd cryptominisat5
$ mkdir build
$ cd build
$ cmake -DUSE_GAUSS=ON -DM4RI_LIBRARIES=../../m4ri/install/lib/libm4ri.a -DM4RI_INCLUDE_DIRS=../../m4ri/install/include -DSTATICCOMPILE=ON ..
$ make
$ cp cryptominisat5 ../../../cryptominisat4
```

## Contacts

* Seonmo Kim (smkim@cs.umn.edu)
* Stephen McCamant (mccamant@cs.umn.edu)

-- INSTALL
# Install STP

STP relies on a few dependencies, namely boost, flex/bison and minisat. Installing the required header files and binaries can be achieved through the following. (Tested on Ubuntu 14.04.)
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

# Install Cryptominisat4 with projection

Make sure to build with Gaussian Eliminiation.
```
$ git clone https://gitlab.com/QIF/cryptominisat4.git
$ cd cryptominisat
$ mkdir build && cd build
$ cmake -DUSE_GAUSS=ON ..
$ make
```
Then, copy the Cryptominisat binary file (``cryptominisat4``) into SearchMC folder.

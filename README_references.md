# References

curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).
```
SPDX-License-Identifier: GPL-3.0-or-later
Copyright (C) 2022-today --- Leonardo de Oliveira Martins
```
curupixa (and other software wrote by me) borrows or was inspired by many existing open source libraries and publicly 
available implementations.
Below follows a reference list of libraries and function implementations which were used, studied, or even considered when writting the software. 
Some minor or well-stablished (public domain) functions or tricks may be missing from here, but are still mentioned as comments in the 
code to the best of my recollection.

## List of libraries and algorithms supporting curupixa

* https://github.com/lemire/testingRNG (Apache License 2.0)
* https://github.com/opencoff/portable-lib (GPL-2.0)
* https://github.com/veorq/SipHash (CC0)
* https://github.com/maciejczyzewski/retter (CC0

## List of other libraries and algorithms 

These libraries are not implemented in curupixa, but have helped me in other software or served as inspiration somehow
(aesthetically, or through OO ideas). The list below is copied from [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib)
which still is my main repository of code.

### hash functions

 * http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200101/homework10/hashfuncs.html
 * http://burtleburtle.net/bob/hash/integer.html 
 * http://www.cse.yorku.ca/~oz/hash.html
 * https://lemire.me/blog/2018/08/15/fast-strongly-universal-64-bit-hashing-everywhere/
 * https://github.com/Cyan4973/smhasher
 * https://github.com/PeterScott/murmur3/ 
 * http://www.concentric.net/~Ttwang/tech/inthash.htm 

The rolling hash algorithm for genomics (not incorporated into biomcmc-lib yet) was inspired by the [linclust algorithm](https://github.com/soedinglab/MMseqs2). 
Their implementation is clean and very fast, but (at least at the time of our implementation) does not compute the
reverse strand &mdash; since it works with a reduced amino acid alphabet. 
 
### hash table

Our first implementation is derived from the software (Rec-I-DCM3)[https://web.njit.edu/~usman/RecIDCM3.html], 
released under the GPL license (Copyright (C) 2004 The University of Texas at Austin. 

### argtable 

The files `argtable3.c` and `argtable3.h` were copied with small modicifcations from [argtable3](https://www.argtable.org/) library, 
maintained by Tom G. Huang at the time of this import (2019.06) and are distributed under a BSD license. 
Please refer to  https://github.com/argtable/argtable3 for the list of authors.

### clustering

 *  OPTICS algorithm based on https://github.com/guineri/GOPTICS  

### quickselect
Code by Nicolas Devillard - 1998. Public domain.  http://ndevilla.free.fr/median

 * Quickselect algorithm: described in "Numerical recipes in C" Section 8.5 (ISBN 0-521-43108-5)
 * Writh algorithm: "Algorithms + data structures = programs". Englewood Cliffs: Prentice-Hall, 1976

### getline

Our implementation is originally from the CvsGui project (http://www.wincvs.org/).
It may be equivalent to the standard GNU getline since (POSIX.1-2008)[https://man7.org/linux/man-pages/man3/getline.3.html].

### Hungarian method

The hungarian method was copied (with modifications) from http://www.informatik.uni-freiburg.de/~stachnis/misc.html
The original message follows:
```
 libhungarian by Cyrill Stachniss, 2004  Solving the Minimum Assignment Problem using the
 Hungarian Method.         ** This file may be freely copied and distributed! **
 Parts of the used code was originally provided by the "Stanford GraphGase", but I made changes to this code.
 As asked by  the copyright node of the "Stanford GraphGase", I hereby proclaim that this file are *NOT* part of the
 "Stanford GraphGase" distribution!
```

### Probability distributions

The probability distribution functions and auxiliary mathematical functions are derived from the (R project for Statistical Computing)[http://www.r-project.org/],
version 2.9.1, available under the GPL license.
The code for the discrete sampling comes from the [GNU Scientific Library](https://www.gnu.org/software/gsl/), version 1.14.

### Random number generation

The basic algorithms for random number generation come from the [GNU Scientific Library](https://www.gnu.org/software/gsl/), version 1.14, with 
motivation for the uncorrelated parallel streams from the [Scalable Parallel Pseudo Random Number Generators Library (SPRNG)](http://www.sprng.org/).
Other references:

* http://www.mcs.anl.gov/~itf/dbpp/text/node119.html (random tree method)
* http://ianbullard.squarespace.com/journal/2009/4/28/why-you-should-never-use-rand.html (gamerand)
* http://www.springerlink.com/content/ek34q1jwla5b6438
* http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps (L'ecuyer list of Tausworthe parameters)
* http://wwwmaths.anu.edu.au/~brent/random.html (generalised Maraglia's xorshift)
* http://prng.di.unimi.it/xoroshiro128plusplus.c and http://xoroshiro.di.unimi.it/xoroshiro128plus.c (xoroshiro128)
* https://github.com/FastFilter/xor_singleheader/blob/master/include/xorfilter.h (splitmix)
* http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html (MT19937)

It is noted that most of these algorithms are in the public domain.

### reconciliation

All reconciliation algorithms were implemented here by Leo Martins, based on methods described in

* https://doi.org/10.1137/S0097539798343362 (SIAM J. Comput, 2000)
* http://dx.doi.org/10.1109/TCBB.2010.14 (IEEE/ACM Trans. Comput. Biol. Bioinform., 2010)
* http://dx.doi.org/10.1371/journal.pcbi.1000501 (PLoS Comput. Biol., 2009)
* http://dx.doi.org/10.1006/mpev.1997.0434 (Mol. Phylogenet. Evol., 1997)
* http://dx.doi.org/10.1093/bioinformatics/17.9.821 (Bioinformatics, 2001)

(this is **not** an exhaustive list)

### suffix tree

Implementation not functional, but these are the projects that I considered and have been studying.
```
git@github.com:JGAntunes/longest-substrings.git
git@github.com:ScottBogen/mcc.git
git@github.com:shysaur/shysaur-suffixtrees.git
git@github.com:mattporritt/suffix_tree.git
git@github.com:Jodh/Ukkonen_Algorithm.git
```
### SVD

Not implemented yet, but I am considering the implementation from [git@github.com:lucasmaystre/svdlibc.git](https://github.com/lucasmaystre/svdlibc).
It contains a fast implementation of SVD matrix decomposition by Doug Rohde.
As of version 1.4, the original files (SVDLIBC) are explicitly available under a BSD license, according to 
[the github repository](https://github.com/lucasmaystre/svdlibc). 

### hyperloglog

Not finished yet, but my implementation is based on code by Ivan Vitjuk https://github.com/ivitjuk/libhll, 
released under an ISC License.

# Source references and attribution

curupixa is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).
```
SPDX-License-Identifier: GPL-3.0-or-later
Copyright (C) 2022-today --- Leonardo de Oliveira Martins
```
The references below record source implementations, algorithm descriptions and historical influences on curupixa. Current module behaviour is described in [Migration coverage](migration.md).

## List of libraries and algorithms supporting curupixa

* https://github.com/lemire/testingRNG (Apache License 2.0)
* http://burtleburtle.net/bob/rand/smallprng.html (CC0)
* http://prng.di.unimi.it/xoroshiro128plusplus.c (CC0)
* https://github.com/opencoff/portable-lib (GPL-2.0)
* https://github.com/veorq/SipHash (CC0)
* https://github.com/maciejczyzewski/retter (CC0)
* https://github.com/mailund/hash/ (GPL-3.0)
* https://graphics.stanford.edu/~seander/bithacks.html (CC0)
* https://prng.di.unimi.it (CC0)
* https://github.com/PeterScott/murmur3 (CC0) 
* https://github.com/wolkykim/qlibc/blob/master/src/utilities/qhash.c (BSD-like)
* https://github.com/yinqiwen/geohash-int/blob/master/geohash.c (BSD3)
* https://github.com/skeeto/hash-prospector (CC0)

* https://www.gnu.org/software/gsl/ (GPL-3.0)
* https://github.com/cran/qrng/ (GPL-3.0)



## Algorithm and historical references
The references below include historical sources from [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) as well as the redesigned implementations. The [migration inventory](migration.md) identifies current capabilities and deliberate changes. Source files retain the attribution for incorporated code; a historical reference does not imply that the current module copies that implementation.

#### hash functions

 * http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200101/homework10/hashfuncs.html
 * http://burtleburtle.net/bob/hash/integer.html 
 * http://www.cse.yorku.ca/~oz/hash.html
 * https://lemire.me/blog/2018/08/15/fast-strongly-universal-64-bit-hashing-everywhere/
 * https://github.com/Cyan4973/smhasher
 * http://www.concentric.net/~Ttwang/tech/inthash.htm (now available at https://gist.github.com/badboy/6267743)
 * http://burtleburtle.net/bob/hash/integer.html

The earlier rolling-hash work was inspired by the [linclust algorithm](https://github.com/soedinglab/MMseqs2). The current API provides a byte rolling hash and a separate canonical DNA k-mer iterator.
Their implementation is clean and very fast, but (at least at the time of my implementation) does not compute the
reverse strand &mdash; since it works with a reduced amino acid alphabet. 
 
#### hash table

The historical biomcmc-lib implementation derives from [Rec-I-DCM3](https://web.njit.edu/~usman/RecIDCM3.html), released under the GPL (Copyright (C) 2004 The University of Texas at Austin). Curupixa now uses independently owned byte keys and open addressing.

#### argtable 

Biomcmc-lib includes the BSD-licensed [argtable3](https://www.argtable.org/) command-line parser. Curupixa leaves this third-party utility outside its phylogenetic API.

#### clustering

 *  OPTICS algorithm based on https://github.com/guineri/GOPTICS  

#### quickselect
Code by Nicolas Devillard - 1998. Public domain.  http://ndevilla.free.fr/median

 * Quickselect algorithm: described in "Numerical recipes in C" Section 8.5 (ISBN 0-521-43108-5)
 * Writh algorithm: "Algorithms + data structures = programs". Englewood Cliffs: Prentice-Hall, 1976

#### getline

The legacy reader credits the CvsGui project (http://www.wincvs.org/). The new stream API uses its own buffered line reader and distinguishes EOF, I/O failure and allocation failure.

#### Hungarian method

The historical Hungarian implementation was copied (with modifications) from http://www.informatik.uni-freiburg.de/~stachnis/misc.html. The current checked rectangular assignment API uses a separate shortest-augmenting-path implementation.
The original message follows:
```
 libhungarian by Cyrill Stachniss, 2004  Solving the Minimum Assignment Problem using the
 Hungarian Method.         ** This file may be freely copied and distributed! **
 Parts of the used code was originally provided by the "Stanford GraphGase", but I made changes to this code.
 As asked by  the copyright node of the "Stanford GraphGase", I hereby proclaim that this file are *NOT* part of the
 "Stanford GraphGase" distribution!
```

#### Probability distributions

The probability distribution functions and auxiliary mathematical functions are derived from the (R project for Statistical Computing)[http://www.r-project.org/],
version 2.9.1, available under the GPL license.
The historical discrete sampler credits the [GNU Scientific Library](https://www.gnu.org/software/gsl/), version 1.14. The current opaque sampler builds a Walker alias table and uses an explicit caller-owned stream.

#### Random number generation

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

#### Split distances and all-root reconciliation
The greedy dSPR and directional reconciliation implementations adapt Leonardo de Oliveira Martins's `topology_splitset.c` and `topology_mrca.c` from [guenomu](https://bitbucket.org/leomrtns/guenomu) and the newer [genefam-dist](https://github.com/leomrtns/genefam-dist), under their GPL-2.0-or-later terms, incorporated here as GPL-3.0-or-later. Curupixa replaces attached scratch fields, the quadratic MRCA cache and all-pairs disagreement bitsets with separate analysis objects. The [distance guide](cladistic-distances.md) records conventions, corrections and heuristic differences.

#### reconciliation

All reconciliation algorithms were implemented here by Leo Martins, based on methods described in

* https://doi.org/10.1137/S0097539798343362 (SIAM J. Comput, 2000)
* http://dx.doi.org/10.1109/TCBB.2010.14 (IEEE/ACM Trans. Comput. Biol. Bioinform., 2010)
* http://dx.doi.org/10.1371/journal.pcbi.1000501 (PLoS Comput. Biol., 2009)
* http://dx.doi.org/10.1006/mpev.1997.0434 (Mol. Phylogenet. Evol., 1997)
* http://dx.doi.org/10.1093/bioinformatics/17.9.821 (Bioinformatics, 2001)

(this is **not** an exhaustive list)

#### suffix tree

The current substring index uses a fresh suffix-array/LCP implementation. The following suffix-tree projects informed earlier investigations:
```
git@github.com:JGAntunes/longest-substrings.git
git@github.com:ScottBogen/mcc.git
git@github.com:shysaur/shysaur-suffixtrees.git
git@github.com:mattporritt/suffix_tree.git
git@github.com:Jodh/Ukkonen_Algorithm.git
```
#### SVD

The legacy sparse implementation comes from [SVDLIBC](https://github.com/lucasmaystre/svdlibc).
It contains a fast implementation of SVD matrix decomposition by Doug Rohde.
As of version 1.4, the original files (SVDLIBC) are explicitly available under a BSD license, according to 
[the github repository](https://github.com/lucasmaystre/svdlibc). 
Curupixa currently uses an independent scaled dense Jacobi backend, including for materialised sparse inputs. It does not retain SVDLIBC's LAS2 implementation or scalability contract.

#### hyperloglog

Earlier work considered Ivan Vitjuk's [libhll](https://github.com/ivitjuk/libhll), released under an ISC licence. The current implementation follows the classic estimator in [Flajolet, Fusy, Gandouet and Meunier (2007)](https://algo.inria.fr/flajolet/Publications/FlFuGaMe07.pdf).

#### File formats and tree construction
The [Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) defines GFF coordinates, attributes and escaping. The parser implements the subset and validation limits described in its public header.
BioNJ follows [Gascuel (1997)](https://academic.oup.com/mbe/article-pdf/14/7/685/11165731/1gasc.pdf); the [tree design](tree-design.md) records corrections to the local legacy implementation.


### Github page
[https://github.com/leomrtns/curupixa](https://github.com/leomrtns/curupixa)

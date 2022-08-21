<img src="recipe/curupixa.png" height="300" align="right" alt="Curupixa">

# Curupixá

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Introduction

This library provides low level functions used by other phylogenetic programs. 
It will be a replacement or auxiliary to [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib),
but currently it is only a testing environment.
Some functions in `biomcmc-lib` are not thread-safe, e.g. the random functions.
Furthermore when used as an API for other languages, it should not exit with error (it should delegate to the calling
program what to do with the exception).

In all likelihood you won't need to install it, since the derived software will have taken care of it &mdash; assuming you
downloaded the derived software from github with `--recursive`.
You are free to add it to your own project, but unfortunately I may be unable to provide adequate support. 
You can [join the discussion](https://github.com/leomrtns/curupixa/discussions), or [report an issue](https://github.com/leomrtns/curupixa/issues) if you
know how to solve it :wink:. 
See [more on that below](##how-to-use-it-in-other-projects).

**Etymology:** Curupixá is the name of a Brazilian tree (_Micropholis venulosa_), also known as salgueiro, Guajará, Grumixava, Grubixa, cubixa, Abiurana, Bacu Mixa, etc. 
## Installation
This library is usually not installed directly, but as a submodule of another project. 
It includes, however, the `makefile.am` and `configure.ac` for autotools, and it provides unit tests from the
[libcheck](https://github.com/libcheck/check) library as well as custom checks. 

For full functionality you will need `openmp`, `libcheck`, `zlib`, `bz2`, and `liblzma`. 
If you can install them system-wide:
```[bash]
## 'bootstrap' the configuration files (needed when cloning from github):
/home/simpson/$ apt-get install pkg-config autotools-dev autoconf automake libtool
/home/simpson/$ (cd curupixa && autoreconf)  ## the parentheses avoid entering the directory afterwards
## install extra libraries: 
/home/simpson/$ apt-get install zlib1g-dev libomp-dev libbz2-dev check liblzma-dev
```

If you do not have admin privileges but you have conda/mamba:
```
/home/simpson/$ mamba install automake libtool pkg-config make libgcc-ng check zlib xz bzip2 libgomp
/home/simpson/$ (cd curupixa && autoreconf)
```

If you really want to install this library separately but you do not have admin privileges and you don't want to use conda, 
then you'll have to pretend we are in the nineties and install all dependencies by hand. 
The libraries rely on `pkg-config` to find their location, otherwise you'll see a cryptic error message about `possibly
undefined macro`. 
If you can only install `pkg-config` through conda then you may need to install the C libraries via conda as well.
Or checking and updating your [`$PKG_CONFIG_PATH` environment variable](https://people.freedesktop.org/~dbn/pkg-config-guide.html).

This sofware uses `autotools`, which means you can install it with `configure` + `make`.
You may need to define where you want it installed with `configure --prefix=DIR`, which is where are your unix-like
`include/`, `lib/`, and `bin/` directories. My favourite is `~/local`. Another popular option if you are on a conda
environment is `--prefix=${CONDA_PREFIX}`

Here is an example of its installation, although **in most cases you do not need to install this library separately**:

```[bash]
/home/simpson/$ git clone https://github.com/leomrtns/curupixa.git
/home/simpson/$ cd curupixa && ./autogen.sh
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../curupixa/configure --prefix=${HOME}/local
/home/simpson/$ make; make install
/home/simpson/$ make check  # battery of unit and integration tests for both tatajuba and biomcmc-lib
```

Notice that two libraries are installed, `libcurupixa.la` and `libcurupixa_static.la` (with `.a` and `.so` versions). 
The first (`libcurupixa.la`) is installed with `make install` and is available "globally" or wherever you tell it with
`--prefix`.
The second version is a [convenience library](https://www.gnu.org/software/automake/manual/html_node/Libtool-Convenience-Libraries.html)
which is used with other programs (through `autotools` subdirs usually) and is not shared. 

Remember that the installation and `autogen.sh` in particular modify/add local files; therefore updating the repository
from github will complain about uncommited changes. You can run `git stash` (or reinstall from scratch) before `git pull`.

## How to use it in other projects 

As I said, this is not a general-purpose library and I can't help much you if you want to use it as such. So in the following
I assume you know enough to be able to include in your own project. If you don't use autotools you can install it globally. 

The main component is to include something like the following in your `configure.ac`
```
AC_SUBST([CURUPIXALIB], [curupixa])
dnl generate makefiles (last recipes, after defining CFLAGS etc.)
AC_CONFIG_FILES([curupixa/Makefile
```
such that the variable `@CURUPIXALIB@` is available to your `makefile.am` as the location of the `curupixa` directories.

There is an older approach, which is to call `curupixa/configure` recursively:
```
AC_CONFIG_SUBDIRS([curupixa])
```
I don't use it anymore due to the issue of autotools+github which forces us to run `autoreconf -fi` (which I provide
through the script `autogen.sh`). And `curupixa` (as well as my other project, `biomcmc-lib`), can generate a
"convenience library" (`noinst_LTLIBRARIES = libcurupixa_static.la`) which is included into the downstreal project
instead of being installed in a shared location.  

Provided that you follow the licence and proper attribution standards, you are free to copy chunks of code and
add into your own projects as an alternative to linking to curupixa. 
Substantial code from external sources are mentioned in the file together with their copyright.

## Acknowledgements 
I incorporated or was inspired by many existing open source libraries and publicly available implementations.
For a more complete, but perhaps non-exhaustive list, please see the file [README_references.md](README_references.md).

## Licence 
SPDX-License-Identifier: [GPL-3.0-or-later](LICENSE)

Copyright (C) 2022-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

The [logo owl](https://www.flickr.com/photos/biodivlibrary/51267563412/in/album-72157719506234015/) is based on the 
[Fauna japonica v.4 Aves](https://www.flickr.com/photos/biodivlibrary/albums/72157719506234015/with/51267563412/)
provided by the Biodiversity Heritage Library under the Public Domain.

![Anurag's github stats](https://github-readme-stats.vercel.app/api?username=leomrtns&count_private=true&show_icons=true&theme=calm)


#!/bin/sh

export AUTOMAKE="automake --foreign -a"

# autogenerated files
rm -f aclocal.m4 m4/lib* m4/lt*
autoreconf -f -i

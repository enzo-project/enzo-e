#!/bin/bash
#
# This script compares Enzo-P's "new" PPM with that in ENZO's "week of code"
# development branch.  ENZOPSRC and ENZOSRC should be set appropriately.
# Output is a set of files, one per Fortran source file, with differences
# between source files (ignoring spacing or newlines)

ENZOPSRC=$HOME/Cello/cello-src/src/Enzo
ENZOSRC=$HOME/Enzo/enzo-dev/src/enzo
FILES="calcdiss calc_eigen euler flux_hll flux_twoshock inteuler intpos intprim intvar pgas2d pgas2d_dual twoshock"

for file in $FILES; do
    echo $file
    diff -bB $ENZOSRC/$file.F $ENZOPSRC/${file}_dev.F > diff.$file
done


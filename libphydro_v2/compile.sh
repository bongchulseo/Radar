#!/bin/sh

LIBNAME="libphydro"
rm $LIBNAME.a
for x in `ls *.c`;
do
  gcc -c $x -m64 -Ofast -flto -march=native -funroll-loops
done

ar rc $LIBNAME.a *.o
ranlib $LIBNAME.a
rm *.o
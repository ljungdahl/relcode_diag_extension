#!/bin/bash

CMAKE=/home/ljungdahl/source/cmake-3.18.1-Linux-x86_64/bin/cmake

cd build

${CMAKE} -D CMAKE_Fortran_COMPILER=gfortran-9 -D CMAKE_C_COMPILER=gcc-9 ..
cd ..

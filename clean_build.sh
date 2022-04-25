#!/bin/bash

CMAKE=/home/ljungdahl/source/cmake-3.18.1-Linux-x86_64/bin/cmake

cd build
rm -rf *
${CMAKE} ..
cd ..

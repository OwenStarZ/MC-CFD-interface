#!/bin/sh
rm -r CMakeFiles
rm CMakeCache.txt
rm Makefile
rm cmake_install.cmake

cmake ..
make

#!/bin/sh
set -e
python init.py
rm -f ./lib_h5rw/lib/*
cd lib_h5rw/build
sh build.sh
cd ../..
cp ./lib_h5rw/lib/libh5rw.so ./
cd libudf
rm -f ./lnamd64/3ddp_host/* ./lnamd64/3ddp_node/*
make
cd ..

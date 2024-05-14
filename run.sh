#!/bin/sh
python gen_udf.py
g++ init.cpp -o init -I/usr/local/hdf5/1.10.5/include -L/usr/local/hdf5/1.10.5/lib -lhdf5
./init
rm ./lib_h5rw/lib/*
cd lib_h5rw/build
sh build.sh
cd ../..
cp ./lib_h5rw/lib/libh5rw.so ./
cd libudf
rm ./lnamd64/3ddp_host/* ./lnamd64/3ddp_node/*
make
cd ..

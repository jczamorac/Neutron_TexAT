#!/bin/sh

cd build
cmake ../
make -j 6
mv ndet_TexAT ../
cd ../

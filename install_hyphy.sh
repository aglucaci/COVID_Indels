#!/bin/bash
#Last update 10-04-2020

#Download and install HyPhy
git clone --branch develop https://github.com/veg/hyphy.git hyphy-develop

cd hyphy-develop

cmake ./
# for silverback
#cmake . -DNOAVX=ON

make -j MP
# for silverback
#make -j MPI

make install

#Download hyphy standalone analyses
cd ..
git clone https://github.com/veg/hyphy-analyses.git

#End of file





#!/bin/bash
#Last update 05-28-2020

#Download and install HyPhy
git clone https://github.com/veg/hyphy.git hyphy-develop

cd hyphy-develop

git checkout 2.5.14
#git checkout origin/develop

cmake ./
make -j MP
make -j MPI

#Download hyphy standalone analyses
cd ..
git clone https://github.com/veg/hyphy-analyses.git

#End of file


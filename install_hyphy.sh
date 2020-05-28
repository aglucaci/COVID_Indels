#!/bin/bash


git clone https://github.com/veg/hyphy.git hyphy-develop



cd hyphy-develop
#git checkout origin/develop



cmake ./
make -j MP
make -j MPI



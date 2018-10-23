#!/usr/bin/env bash

conda env create -f `dirname "$0"`/environment.yml
source activate rcd-env
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/x2z .
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make install
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/MESS .
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make install
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/EStokTP .
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc 
make install
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/QTC .
pip install .

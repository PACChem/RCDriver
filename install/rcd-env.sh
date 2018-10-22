create -f `dirname "$0"`/envs/enironment.yml
source activate rcd-env
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/x2z .
cmake -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
make install
cd $(mktemp -d)
git clone --recursive https://github.com/PACChem/QTC .
pip install .

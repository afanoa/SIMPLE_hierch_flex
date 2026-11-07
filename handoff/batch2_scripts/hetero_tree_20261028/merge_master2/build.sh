#!/bin/bash
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/build
hostname > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/build_host.txt; date > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/build_started.txt
cmake /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/src -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DUSE_ARCHOPT=ON -DUSE_OPENMP=ON > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/cmake.log 2>&1 && make -j 32 install > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/make.log 2>&1
echo $? > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/build.exit; date > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/build_finished.txt

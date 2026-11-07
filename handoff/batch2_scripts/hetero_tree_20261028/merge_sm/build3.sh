#!/bin/bash
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_sm
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:/mnt/nasapps/production/python/3.12/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd $B/build && make -j 32 install > $B/make3.log 2>&1; echo $? > $B/build3.exit; date > $B/build3_finished.txt

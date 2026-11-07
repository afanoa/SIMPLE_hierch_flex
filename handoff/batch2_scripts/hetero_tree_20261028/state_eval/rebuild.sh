#!/bin/bash
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
rm -f /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/build2.exit /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/fast_gate2.exit
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/src && git fetch -q ~/SIMPLE state-eval-in && git checkout -q --detach FETCH_HEAD || { echo 9 > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/build2.exit; exit 1; }
git log --oneline -1 | cut -c1-90 > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/rebuild_commit.txt
cd /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/build && make -j 32 install > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/make2.log 2>&1; echo $? > /mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/build2.exit

#!/bin/bash
# rebuild merge_master2 at the scratch ref's tip (incremental), then submit the fast gate on its build tree
D2=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd $D2/src && git fetch -q ~/SIMPLE tree-master-merge-in && git checkout -q --detach FETCH_HEAD && git log --oneline -1 > $D2/rebuild_commit.txt
cd $D2/build && make -j 32 install > $D2/make2.log 2>&1; echo $? > $D2/build2.exit
[ "$(cat $D2/build2.exit)" = 0 ] || exit 1
cat > $D2/fast_gate.sbatch <<EOS
#!/bin/bash
#SBATCH -J fastgate_mm
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH --output=$D2/fast_gate.out
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/nasapps/production/cmake/4.0.3/bin:/mnt/nasapps/production/python/3.12/bin:\$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
cd $D2/src && scripts/run_fast_gate.sh $D2/build > $D2/fast_gate.log 2>&1; echo \$? > $D2/fast_gate.exit
EOS
sbatch --parsable $D2/fast_gate.sbatch > $D2/fast_gate.jobid

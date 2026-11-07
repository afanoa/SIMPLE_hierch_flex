#!/bin/bash
#SBATCH -J seltest
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:10:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/seltest_%A.out
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/merge_master2/install/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
D=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/seltest; mkdir -p $D; cd $D
for k in qtoy2:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates/mm2_qtoy2/1_hetero_tree/node_r qtoynull:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates/mm2_qtoynull/1_hetero_tree/node_r olenull:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates/mm2_olenull/1_hetero_tree/node_r r131root:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates/mm2_r131/1_hetero_tree/node_r r131:/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/gates/mm2_r131/1_hetero_tree/node_r1; do
  n=${k%%:*}; N=${k#*:}; cp $N/n.simple $n.simple
  simple_exec prg=selection projfile=$n.simple oritype=ptcl3D infile=$N/seeds/seeds.txt mkdir=no > $n.log 2>&1; echo "$n selection exit $? $(grep -c 'Inconstitent' $n.log)"
  simple_exec prg=prune_project projfile=$n.simple mkdir=no > ${n}_prune.log 2>&1; echo "$n prune_project exit $? $(grep -ch 'ERROR' ${n}_prune.log)"
done

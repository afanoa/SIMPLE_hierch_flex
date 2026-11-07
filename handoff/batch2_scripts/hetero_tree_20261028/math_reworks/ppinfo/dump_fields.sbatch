#!/bin/bash
#SBATCH -J ppinfo_dump
#SBATCH -p csbdevel
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/logs/%x_%j.out
# dump the stk and ptcl3D segments of a node project and its stage-B project to text: NODE=<node dir> TAG=<name>
set -uo pipefail
: ${NODE:?} ${TAG:?}
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/detlim2/build_vd/install
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$B/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 SIMPLE_PATH=$B OMP_NUM_THREADS=1
D=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/ppinfo/dumps/$TAG; mkdir -p $D; cd $D
for seg in stk ptcl3D ptcl2D; do simple_exec prg=print_project_field projfile=$NODE/n.simple oritype=$seg mkdir=no > n_$seg.txt 2> n_$seg.err; done
simple_exec prg=print_project_field projfile=$NODE/stageB/b.simple oritype=ptcl3D mkdir=no > b_ptcl3D.txt 2> b_ptcl3D.err
simple_exec prg=print_project_info projfile=$NODE/n.simple mkdir=no > n_info.txt 2>&1
wc -l *.txt; head -2 n_ptcl3D.txt; head -2 b_ptcl3D.txt; head -2 n_stk.txt

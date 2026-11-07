#!/bin/bash
#SBATCH -p norm
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:15:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# Which per-particle / per-micrograph quality fields do the RNase P projects carry (a corr-free quality covariate)?
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
P=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/gate/m33_root_m79/1_hetero_tree/node_r/recB0/b.simple
for o in ptcl3D ptcl2D stk mic; do echo "== $o"; simple_exec prg=print_project_field projfile=$P oritype=$o 2>/dev/null | grep "=" | head -2 | cut -c1-900; done
R=/mnt/beegfs/elmlund/afan/rnasep_klem/stree_port/full692_design80_s1/rnasep.simple; [ -s $R ] && for o in mic stk; do echo "== $R $o"; simple_exec prg=print_project_field projfile=$R oritype=$o 2>/dev/null | grep "=" | head -2 | cut -c1-900; done

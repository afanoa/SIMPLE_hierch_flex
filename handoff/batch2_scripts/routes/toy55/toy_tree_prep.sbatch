#!/bin/bash
#SBATCH -J toyah_tree_prep
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Run the real-data chain on the apo/holo toy: meta npz from the toy's pose table, a project copy with the compute
# environment set, then the chain (RTAG toyah_tree) with modest parallelism, then the truth-coloured global embedding.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/toyah_tree; mkdir -p $R; cd $R
$PY $W/make_toy_meta.py $W/lat_ah/ptcl3D_r0.txt $R/meta_toy.npz
cp -f $W/projah/1_fastposes/projah.simple $R/src.simple
simple_exec prg=update_project projfile=$R/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=5 nthr=8 mkdir=no > up.log 2>&1 && echo "project copy with compute environment: $R/src.simple"
cd $W && sbatch --export=ALL,SMOKE=0,RTAG=toyah_tree,SRC=$R/src.simple,METAALL=$R/meta_toy.npz,NP=5,NPARTS=4,NTHR=8,PARMAX=1 -J toyah_tree f692_chain6.sbatch

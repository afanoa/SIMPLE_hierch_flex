#!/bin/bash
#SBATCH -J synth_tree_prep
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Dev chain (chain9, DEV=1, board server) on a synthetic benchmark (env MODE arms|tethered, from run_toysynth_sim):
# meta npz from the simulation's refined pose table, project copy with the compute environment, then the chain.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; R=$A/rnasep_klem/synth_$MODE; mkdir -p $R; cd $R
$PY $W/make_toy_meta.py $W/lat_$MODE/ptcl3D_r0.txt $R/meta_toy.npz || exit 1
RD=$(ls -td $W/proj$MODE/*_refine3D | head -1); cp -f $RD/proj$MODE.simple $R/src.simple
simple_exec prg=update_project projfile=$R/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > up.log 2>&1 && echo "project copy with compute environment: $R/src.simple (from $RD)"
cd $W && sbatch -c 48 --mem=120G --export=ALL,RTAG=synth_$MODE,DEV=1,SRC=$R/src.simple,METAALL=$R/meta_toy.npz -J synth_${MODE}_dev f692_chain9.sbatch

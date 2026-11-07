#!/bin/bash
#SBATCH -J f692sw_prep
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the swing synthetic through the f692e route (baseline organization): waits for the simulation job SIMJ, then meta npz from the
# near-true pose table, a project copy for distributed workers, and chain9 (routes/toy55_f692v2 = f692e route + guards, csbdevel)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; W=$A/routes/toy55_f692v2; TAG=${TAG:-swing30}; RT=${RTAG:-f692sw30}; R=$A/rnasep_klem/$RT
[ -n "${SIMJ:-}" ] && while squeue -h -j $SIMJ 2>/dev/null | grep -q .; do sleep 60; done
[ -s $T/lat_$TAG/ptcl3D_r0.txt ] && [ -s $T/proj$TAG/proj$TAG.simple ] || { echo "simulation output missing"; exit 1; }
mkdir -p $R; cd $R; $PY $T/make_toy_meta.py $T/lat_$TAG/ptcl3D_r0.txt $R/meta_toy.npz || exit 1
cp -f $T/proj$TAG/proj$TAG.simple $R/src.simple; simple_exec prg=update_project projfile=$R/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > up.log 2>&1 || { echo "update_project FAILED"; exit 1; }
cd $W && sbatch -p csbdevel -c 8 -J ${RT}_chain --export=ALL,RTAG=$RT,SRC=$R/src.simple,METAALL=$R/meta_toy.npz,DEV=${DEV:-1},LOCALDEV=0,SERVER=1,NPARTS=10,NP=10,PARMAX=6,CTFABS=0,PSCALE=0,REFINE_EXTRA=ml_reg=no,NFREE_SUB=1,SIMPLE_QSYS_PARTITION=csbdevel f692_chain9.sbatch

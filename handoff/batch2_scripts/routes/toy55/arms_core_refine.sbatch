#!/bin/bash
#SBATCH -J synth_corepose
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# re-refine the poses of a hinged-arm synthetic set against the CORE only (spherical mask of the core diameter, coarse band)
# starting from the chain's refined project, then the pose error against the truth. env MODE (arms3), MSK (122), LP (20), MAXITS (5), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:-arms3}; MSK=${MSK:-122}; NP=${NP:-10}
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
O=$W/corepose_$MODE${TAG:-}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') corepose $MODE] $*"; }
[ -s p.simple ] || cp $A/rnasep_klem/synth_$MODE/src.simple p.simple
simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=14000 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
log "refine3D on the core only: mskdiam $MSK, lp ${LP:-20}, maxits ${MAXITS:-5}, $NP parts x 8"; t0=$(date +%s)
simple_exec prg=refine3D projfile=p.simple vol1=${VOL:-$W/maps160/synth_arms/consensus_synth.mrc} nstates=1 refine=${REFINE:-prob} objfun=cc pgrp=c1 mskdiam=$MSK lp=${LP:-20} maxits=${MAXITS:-5} nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-160)"; exit 1; }
RD=$(ls -td $O/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $O/ptcl3D_core.txt
log "refined in $(( $(date +%s) - t0 )) s -> $RD/p.simple, $(grep -c . $O/ptcl3D_core.txt) rows"
$PY $W/pose_err_hist.py $W/sim$MODE/all.txt $O/ptcl3D_core.txt $W/lat_$MODE/ptcl3D_r0.txt $W/sim$MODE/gt_particles.txt 2>&1 | grep -v Warning | sed 's/^start /chain-refined/'
log DONE

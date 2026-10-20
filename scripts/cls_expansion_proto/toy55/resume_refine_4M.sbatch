#!/bin/bash
#SBATCH -J toy55_scale
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH -t 48:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/toy55_scale.%j.out
# resume the 4M set at the refine3D step (the import is done; the 48 GB master was OOM-killed at 130 GB RSS)
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=55_4M; NP=${NP:-10}; PJ=proj$MODE; S=$W/sim$MODE; log(){ echo "[$(date '+%F %T') scale $MODE] $*"; }; cd $W/$PJ
simple_exec prg=update_project projfile=$PJ.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=86000 nparts=$NP nthr=8 mkdir=no > up2.log 2>&1
log "refine3D from the toy consensus, perturbed start (angerr 10), $NP parts (resumed, 200 GB master)"
simple_exec prg=refine3D projfile=$PJ.simple vol1=$W/maps160/toy_consensus.mrc nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=200 lp=8 maxits=4 nparts=$NP nthr=8 > refine2.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 'ERROR\|THROW' refine2.log | cut -c1-160)"; exit 1; }
RD=$(ls -td $W/$PJ/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $W/lat_$MODE/ptcl3D_r0.txt
log "DONE synthetic $MODE: $S, project $RD/$PJ.simple, pose table $W/lat_$MODE/ptcl3D_r0.txt ($(grep -c . $W/lat_$MODE/ptcl3D_r0.txt) rows)"

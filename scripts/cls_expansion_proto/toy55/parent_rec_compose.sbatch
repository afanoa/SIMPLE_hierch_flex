#!/bin/bash
#SBATCH -J klem_parent_compose
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the parent for core borrowing = an UNFILTERED even/odd reconstruction of all refined particles at the consensus poses
# (the deposited-style rec_final halves correlate to Nyquist, so they are not independent half-sets); then the composed
# state maps of the production leaf set and the zoomed sub-states with FSC readouts (region weighted by mask^2)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; O=$A/rnasep_klem/parent692; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') parent] $*"; }
if [ ! -s $O/parent_even.mrc ]; then
  cp -f $A/rnasep_coupled/full692k_2d/rnasep.simple p.simple; $PY -c "import numpy as np; fl=np.loadtxt('$A/rnasep_coupled/flags_refined.txt',dtype=int); np.savetxt('$O/states.txt', fl, fmt='%d'); print('refined rows', fl.sum())"
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=20000 nparts=10 nthr=8 mkdir=no > up.log 2>&1
  log "reconstruct3D of all refined particles at the consensus poses (10 parts x 8)"; simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=10 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
  RR=$(ls -td $O/*_reconstruct3D | head -1); cp -f $RR/recvol_state01_even.mrc $O/parent_even.mrc; cp -f $RR/recvol_state01_odd.mrc $O/parent_odd.mrc; cp -f $RR/recvol_state01.mrc $O/parent.mrc; rm -rf $RR/*_part*.mrc 2>/dev/null
fi
for C in consolidated_12_fullw_region consolidated_12_fullw_zoom; do
  echo "== f692e/$C composed with the 692k parent (unfiltered halves)"; $PY $W/compose_state_maps.py $A/rnasep_klem/f692e $A/rnasep_klem/f692e/$C $O/parent $A/rnasep_klem/f692e/$C/composed 2>&1 | grep -v Warning
done
echo "COMPOSE2 DONE"

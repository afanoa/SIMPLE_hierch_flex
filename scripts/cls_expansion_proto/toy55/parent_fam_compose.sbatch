#!/bin/bash
#SBATCH -J klem_parentfam
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-family parents: unfiltered even/odd reconstruction of ALL particles of family F at that family's refined poses
# (the states of a family are reconstructed in that frame; the consensus-frame parent mismatched family 2), then the
# composed state maps (family core outside the state's region, the state inside) with FSC readouts
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692e; O=$R/parents; mkdir -p $O; log(){ echo "[$(date '+%F %T') parentfam] $*"; }
for F in 1 2; do
  [ -s $O/fam${F}_even.mrc ] && continue; D=$O/rec$F; mkdir -p $D; cd $D; PS=$(ls -td $R/fam$F/refine/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple
  simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s=0; for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); s=a[2]}; print (s>0)?1:0}' > $D/states.txt; log "family $F: $(grep -c 1 $D/states.txt) particles with a state"
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=20000 nparts=10 nthr=8 mkdir=no > up.log 2>&1
  log "reconstruct3D family $F at its refined poses"; simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=1 nparts=10 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
  RR=$(ls -td $D/*_reconstruct3D | head -1); cp -f $RR/recvol_state01_even.mrc $O/fam${F}_even.mrc; cp -f $RR/recvol_state01_odd.mrc $O/fam${F}_odd.mrc; cp -f $RR/recvol_state01.mrc $O/fam$F.mrc; rm -rf $RR/*_part*.mrc 2>/dev/null
done
for C in consolidated_12_fullw_region consolidated_12_fullw_zoom; do
  echo "== f692e/$C composed with the per-family parents"; $PY $W/compose_state_maps.py $R $R/$C "$O/fam{F}" $R/$C/composed 2>&1 | grep -v Warning
done
echo "COMPOSE3 DONE"

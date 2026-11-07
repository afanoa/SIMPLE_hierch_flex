#!/bin/bash
#SBATCH -J igm2_truepose
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The truth-free factorised route on an IgM-like set at NEAR-TRUE poses (the simulation's poses perturbed by ANGERR deg) in
# ONE job: project from the perturbed truth (eo by row parity, state 1) -> consensus reconstructed at those poses -> the
# factorised stages (OSFX=_tp<ANGERR>, PSRC, CONSMAP; junk from the existing tree root) -> readout next to the refined-pose
# and core-realigned runs. Separates the pose loss from the per-Fab signal limit. env MODE (armsigm2), ANGERR (3), SHERR (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KL=$A/rnasep_klem; export MODE=${MODE:-armsigm2}; AE=${ANGERR:-3}
S=$W/sim$MODE; P=$W/projtp_${MODE}_$AE; mkdir -p $P; log(){ echo "[$(date '+%F %T') igm2_truepose $MODE] $*"; }
step(){ local s=$1; shift; local L=$W/igm2_${MODE}_tp${AE}_${s%.sbatch}.log; log "stage $s"; env "$@" bash $W/$s > $L 2>&1; local e=$?; grep -a -E "junk removed|truth check|regions|merged region|gap threshold|FAILED" $L | tail -14; log "stage $s exit $e"; return $e; }
# ---- project at the perturbed true poses
if [ ! -s $P/proj.simple ]; then
  cd $P && simple_exec prg=new_project projname=proj > np.log 2>&1; cd $P/proj || exit 1
  simple_exec prg=oriops oritab=$S/all.txt outfile=$P/all_tp.txt angerr=$AE sherr=${SHERR:-1} mkdir=no > $P/oriops.log 2>&1 || { log "oriops FAILED"; exit 1; }
  awk '{ sub(/ state=[^ ]*/, ""); if ($0 ~ / eo=/) sub(/ eo=[^ ]*/, " eo=" ((NR - 1) % 2)); else $0 = $0 " eo=" ((NR - 1) % 2); print $0 " state=1" }' $P/all_tp.txt > $P/all_tp_import.txt
  simple_exec prg=import_particles projfile=proj.simple stk=$S/all.mrcs oritab=$P/all_tp_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > $P/import.log 2>&1 || { log "import FAILED: $(grep -a -m2 ERROR $P/import.log | cut -c1-160)"; exit 1; }
  mv proj.simple $P/proj.simple
  simple_exec prg=update_project projfile=$P/proj.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=43200 nparts=10 nthr=8 mkdir=no > $P/up.log 2>&1
fi
simple_exec prg=print_project_field projfile=$P/proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $P/ptcl3D_tp.txt
echo "== near-true poses: $($PY $W/pose_error.py $P/ptcl3D_tp.txt $S/all.txt --tag "truth + $AE deg" 2>&1 | tail -1); odd-half rows $(grep -c ' eo=1\.' $P/ptcl3D_tp.txt) of $(grep -c . $P/ptcl3D_tp.txt)"
# ---- consensus at those poses
cd $P; if [ ! -s $P/cons_tp.mrc ]; then cp -f $P/proj.simple $P/p_rec.simple
  simple_exec prg=reconstruct3D projfile=p_rec.simple pgrp=c1 mskdiam=232 nparts=10 nthr=8 > rec.log 2>&1 && RR=$(ls -td $P/*_reconstruct3D | head -1) && cp -f $RR/recvol_state01.mrc $P/cons_tp.mrc || log "reconstruct3D FAILED: $(grep -a -m2 'ERROR\|THROW' rec.log | cut -c1-160)"; fi
CONS=$P/cons_tp.mrc; [ -s $CONS ] || { CONS=$W/proj$MODE/1_refine3D/recvol_state01.mrc; log "fallback consensus: the refined-pose map"; }
RP=$KL/synth_${MODE}_tp$AE; mkdir -p $RP
[ -s $RP/src.simple ] || { cp -f $P/proj.simple $RP/src.simple; simple_exec prg=update_project projfile=$RP/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > $RP/up.log 2>&1; }
# ---- the factorised route at the near-true poses
X="OSFX=_tp$AE PSRC=$RP/src.simple CONSMAP=$CONS"
step igm_factorised.sbatch $X V1BOARDS=0 || exit 1
step igm_factorised2.sbatch $X || exit 1
step igm_factorised3.sbatch $X SKIPMERGE=1 || exit 1
step igm_factorised4.sbatch $X || exit 1
echo "== final readout, near-true poses (truth + $AE deg)"; sed -n '/start selection/,$p' $W/igm2_${MODE}_tp${AE}_igm_factorised4.log | grep -v "^\["
echo "== comparison, refined poses"; grep -a "per arm, best region" $W/igm2_${MODE}_igm_factorised4.log | tail -1
echo "== comparison, core-realigned poses"; grep -a "per arm, best region" $W/igm2_${MODE}_ra_igm_factorised4.log | tail -1
log DONE

#!/bin/bash
#SBATCH -J igm2_realign
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Pose diagnosis and rigid-core realignment of an IgM-like set, then the truth-free factorised route on the realigned poses,
# in ONE job (truth only in readouts): (1) per-image ceiling (arm_llr) at the TRUE poses (lp 20, 10) and the refined poses
# (lp 20); (2) rigid-core mask from the consensus (periphery_mask --core_out); (3) refine3D with pcg_mskfile = that core
# (global prob search, lp 8) on a fresh copy of the untouched full-set refinement project -> pose error before / after;
# (4) ceiling at the realigned poses; (5) full consensus at the realigned poses (reconstruct3D); (6) the factorised stages
# with OSFX=_ra, PSRC = realigned project, CONSMAP = realigned consensus; final readout. env MODE (armsigm2), CORE (0.40), RMAX (5)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KL=$A/rnasep_klem; export MODE=${MODE:-armsigm2}
S=$W/sim$MODE; MK=$W/maps160/synth_$MODE/masks; RA=$W/realign_$MODE; mkdir -p $RA; cd $RA
CONS0=$W/proj$MODE/1_refine3D/recvol_state01.mrc; P0=$W/proj$MODE/1_refine3D/proj$MODE.simple; log(){ echo "[$(date '+%F %T') igm2_realign $MODE] $*"; }
ceil(){ local tag=$1 lp=$2 oris=${3:-}; echo "== ceiling at the $tag poses, lp $lp"; OMP_NUM_THREADS=48 OPENBLAS_NUM_THREADS=48 $PY $W/arm_llr.py $S $W/truth_$MODE.txt $MK/all_arms.mrc $W/maps160/synth_$MODE $W/armllr_$MODE/${tag}_lp$lp --nsub 3000 --lp $lp --far 30 ${oris:+--oris $oris} 2>&1 | grep -v Warning | grep -a -E "accuracy|d' median|arm +[0-9]+:" | head -20; }
step(){ local s=$1; shift; local L=$W/igm2_${MODE}_ra_${s%.sbatch}.log; log "stage $s"; env "$@" bash $W/$s > $L 2>&1; local e=$?; grep -a -E "junk removed|truth check|regions|merged region|gap threshold|FAILED" $L | tail -14; log "stage $s exit $e"; return $e; }
# ---- (1) ceiling at true and refined poses
ceil true 20; ceil true 10; ceil refined 20 $W/lat_$MODE/ptcl3D_r0.txt
# ---- (2) rigid-core mask (data-derived)
$PY $W/periphery_mask.py $CONS0 $RA/periphery0.mrc --core ${CORE:-0.40} --core_out $RA/core.mrc --check $MK 2>&1 | grep -v Warning
# ---- (3) core-masked global refinement on a fresh copy of the untouched refinement project
if [ ! -s $RA/ptcl3D_ra.txt ]; then
  cp -f $P0 $RA/p.simple; simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=43200 nparts=10 nthr=8 mkdir=no > up.log 2>&1
  log "refine3D: pcg_mskfile = data-derived core, global prob search, lp 8, maxits ${RMAX:-5}"
  simple_exec prg=refine3D projfile=p.simple vol1=$CONS0 pcg_mskfile=$RA/core.mrc nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=220 lp=8 maxits=${RMAX:-5} nparts=10 nthr=8 > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-160)"; exit 1; }
  [ "$(grep -a -c 'RESOLUTION @ FSC=0.143' refine.log)" -ge 2 ] || { log "refine3D did not iterate: $(grep -a -m2 'ERROR' refine.log | cut -c1-160)"; exit 1; }
  RD=$(ls -td $RA/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $RA/ptcl3D_ra.txt
fi
RD=$(ls -td $RA/*_refine3D | head -1)
echo "== pose error before / after the core realignment"; $PY $W/pose_error.py $W/lat_$MODE/ptcl3D_r0.txt $S/all.txt --tag "before (full reference)" 2>&1 | tail -1; $PY $W/pose_error.py $RA/ptcl3D_ra.txt $S/all.txt --tag "after (core, pcg_mskfile)" 2>&1 | tail -1
$PY $W/pose_err_hist.py $S/all.txt $RA/ptcl3D_ra.txt $W/lat_$MODE/ptcl3D_r0.txt $S/gt_particles.txt 2>&1 | grep -v Warning | sed 's/^refined /after  /; s/^start /before /'
# ---- (4) ceiling at the realigned poses
ceil realigned 20 $RA/ptcl3D_ra.txt
# ---- (5) full consensus at the realigned poses, realigned project for the boards
if [ ! -s $RA/cons_ra.mrc ]; then cp -f $RD/p.simple $RA/p_ra.simple
  simple_exec prg=reconstruct3D projfile=p_ra.simple pgrp=c1 mskdiam=232 nparts=10 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
  RR=$(ls -td $RA/*_reconstruct3D | head -1); cp -f $RR/recvol_state01.mrc $RA/cons_ra.mrc; fi
RP=$KL/synth_${MODE}_ra; mkdir -p $RP; [ -s $RP/src.simple ] || { cp -f $RD/p.simple $RP/src.simple; simple_exec prg=update_project projfile=$RP/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > $RP/up.log 2>&1; }
# ---- (6) the factorised route on the realigned poses
X="OSFX=_ra PSRC=$RP/src.simple CONSMAP=$RA/cons_ra.mrc"
step igm_factorised.sbatch $X V1BOARDS=0 || exit 1
step igm_factorised2.sbatch $X || exit 1
step igm_factorised3.sbatch $X SKIPMERGE=1 || exit 1
step igm_factorised4.sbatch $X || exit 1
echo "== final readout, realigned poses"; sed -n '/merged stage readout/,$p' $W/igm2_${MODE}_ra_igm_factorised4.log | grep -v "^\["
echo "== for comparison, refined poses (igm2_pipeline)"; grep -a "per arm, best region" $W/igm2_${MODE}_igm_factorised4.log | tail -1
log DONE

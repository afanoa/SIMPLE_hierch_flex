#!/bin/bash
#SBATCH -J igm2_pipeline
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=160G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# A harder IgM-like set through the truth-free factorised route, unchanged, in ONE job (truth only in readouts):
# run_igm_sim (Fab pairs at real IgM ranges, simulation, full-set refinement; NOCHAIN) -> pair-aware sector masks
# (readouts) -> igm_dtree (independent halves, tree ROOT only: the route uses its junk split) -> igm_factorised (junk +
# coarse board; V1BOARDS 0) -> igm_factorised2 (periphery, 14 -> 8 A board, regions) -> igm_factorised3 (residual stage;
# SKIPMERGE) -> igm_factorised4 (merge gap:1.5, 4 starts, selection, readout). Each stage logs to igm2_<MODE>_<stage>.log.
# env MODE (armsigm2), IGM_EXTRA (--pairs --splay 18 --theta 40 --phi 75)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; export MODE=${MODE:-armsigm2}
log(){ echo "[$(date '+%F %T') igm2 $MODE] $*"; }
step(){ local s=$1; shift; local L=$W/igm2_${MODE}_${s%.sbatch}.log; log "stage $s"; env "$@" bash $W/$s > $L 2>&1; local e=$?
  grep -a -v "^ *$" $L | grep -a -E "==|readout|truth|R\^2|corr|region|junk|pose error|refined vs|merged|gap|FAILED|DONE|exit" | tail -60; log "stage $s exit $e"; return $e; }
step run_igm_sim.sbatch NOCHAIN=1 "IGM_EXTRA=${IGM_EXTRA:---pairs --splay 18 --theta 40 --phi 75}" || exit 1
$PY $W/igm_arm_masks.py $W/maps160/synth_$MODE/masks --pairs 2>&1 | grep -v Warning
step igm_dtree.sbatch MAXDEPTH=0 TAG=dtree || exit 1
step igm_factorised.sbatch V1BOARDS=0 || exit 1
step igm_factorised2.sbatch || exit 1
step igm_factorised3.sbatch SKIPMERGE=1 || exit 1
step igm_factorised4.sbatch || exit 1
echo "== final readout (igm_factorised4)"; sed -n '/merged stage readout/,$p' $W/igm2_${MODE}_igm_factorised4.log
log DONE

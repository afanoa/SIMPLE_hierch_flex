#!/bin/bash
#SBATCH -J igm2_avgs
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Sub-class-level E-step in the per-region boards of an IgM-like set (truth only in readouts): the full expansion's 1940
# averages as the E-step units (igm2_avgboard.sbatch array over the merged regions): near-true poses with avg modes pool,
# unmix, unmixc (parent-centred), refined poses with pool; factorised_readout.py per run next to the per-particle boards
# (igm2_bags.sbatch run_tp3_beta0 = same Python board, no pooling). env MODE (armsigm2), CONC (12)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; export MODE=${MODE:-armsigm2}
B=$W/bags_$MODE; O=$W/avgs_$MODE; mkdir -p $O; AVGS=$W/route2d_$MODE/expansion_all/cls_expansion_cavgs.mrc; log(){ echo "[$(date '+%F %T') igm2_avgs $MODE] $*"; }
[ -s $B/ptcl2D_all.txt ] && [ -s $B/meta_tp3.npz ] && [ -s $AVGS ] || { log "needs igm2_bags.sbatch outputs ($B) and $AVGS"; exit 1; }
CFG=$O/tasks.txt; : > $CFG
for spec in "tp3:$W/factorised3_${MODE}_tp3/merged4/regions:$W/projtp_${MODE}_3/cons_tp.mrc:pool unmix unmixc" "ref:$W/factorised3_$MODE/merged4/regions:$W/proj$MODE/1_refine3D/recvol_state01.mrc:pool"; do
  IFS=: read tag rdir cons modes <<< "$spec"
  for m in $modes; do R=$O/run_${tag}_$m; mkdir -p $R/regions; for f in $(ls $rdir/region_*.mrc); do r=$(basename $f .mrc); r=${r#region_}; ln -sfn $f $R/regions/region_$r.mrc; mkdir -p $R/arm_$r
    echo "$R/arm_$r/board $f $B/meta_$tag.npz $cons $m $B/ptcl2D_all.txt $B/parents.txt $AVGS" >> $CFG; done; done; done
N=$(grep -c . $CFG); log "$N boards (array, ${CONC:-12} at a time)"
J=$(sbatch --parsable --array=1-$N%${CONC:-12} --export=ALL,CFG=$CFG,MODE=$MODE $W/igm2_avgboard.sbatch); while squeue -h -j $J 2>/dev/null | grep -q .; do sleep 30; done
log "boards done: $(cat $W/igm2_avgboard.${J}_*.out 2>/dev/null | grep -c 'task done') of $N; failed: $(cat $W/igm2_avgboard.${J}_*.out 2>/dev/null | grep -c FAILED)"
echo "== class-frame calibration and sub-class E-step (first tasks)"; for t in 1 2; do cat $W/igm2_avgboard.${J}_$t.out 2>/dev/null | grep -v "task done"; done
for R in $O/run_tp3_pool $O/run_tp3_unmix $O/run_tp3_unmixc $O/run_ref_pool; do echo "== $(basename $R)"; $PY $W/factorised_readout.py $R $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning | grep -a "per arm"; done
echo "== per-particle Python boards, near-true poses (no pooling)"; $PY $W/factorised_readout.py $B/run_tp3_beta0 $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $W/maps160/synth_$MODE/masks 2>&1 | grep -v Warning | grep -a "per arm"
log DONE

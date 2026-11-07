#!/bin/bash
#SBATCH -J klem_calE_sw2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# More switches for the switch / swing call: cdisc closed along the phase-B score with junk out of the bins (tree calls; truth junk
# as an oracle), and cdisc B (the blob) along its best population-wide score (pick_score_any.py) with tree junk out; 5 at-node bins,
# the root's recB0 reconstruct3D with nstates=5, postprocess=no; phaseE_maptests.py on each
set -uo pipefail
NP=2; NT=4; TAG=calE_sw2; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/switch2; mkdir -p $OUT; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; PB=$K/factorial/phaseB; ST=$K/factorial/junk_test/simsites_truth_jt.txt
for ex in tree truth; do newdir cdisc_closed_${ex}junk; cp $PB/cdisc_closed/work/score.txt .
  { $PY $T/motion_bins4.py score.txt states.txt cdisc $ST $ex $CD > ../bins.txt 2>&1 \
    && recon cdisc_closed_${ex}junk $CD/node_r/stageB/b.simple $CD/node_r/recB0/hetero_step.sh $CD/node_r/mask.mrc closed; } || log "cdisc_closed_${ex}junk FAILED: $(tail -1 ../bins.txt)"
done
newdir cdisc_B_treejunk
{ $PY $T/pick_score_any.py score.txt B_blob $ST $CD > ../pick.txt 2>&1 && $PY $T/motion_bins5.py score.txt states.txt tree $CD $ST A_domain B_blob C_closed > ../bins.txt 2>&1 \
  && recon cdisc_B_treejunk $CD/node_r/stageB/b.simple $CD/node_r/recB0/hetero_step.sh $CD/node_r/mask.mrc B_blob; } || log "cdisc_B_treejunk FAILED: $(tail -qn1 ../pick.txt ../bins.txt 2>&1)"
finish

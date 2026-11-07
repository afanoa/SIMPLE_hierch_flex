#!/bin/bash
#SBATCH -J klem_calE_sw
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Switch calibration of the switch / swing call: CAHRA open along a brightness-corrected score (pick_score_bc.py) and cdisc A with
# junk out of the bins (the tree's own junk calls; the truth junk as an oracle), 5 at-node bins (motion_bins4.py), the root's recB0
# reconstruct3D with nstates=5, postprocess=no (seriesE_lib.sh); phaseE_maptests.py on each and on the three phase-B switch series
set -uo pipefail
NP=2; NT=4; TAG=calE_sw; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/switch; mkdir -p $OUT/ex_switches; CA=$S/arm_M_c11_s1/1_hetero_tree; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; PB=$K/factorial/phaseB
GT=/mnt/beegfs/elmlund/EMPIAR_downloads/13991/data/source.star; ST=$K/factorial/junk_test/simsites_truth_jt.txt
( $PY $T/phaseE_maptests.py ex_cahra_open:$PB/cahra_open:$CA/node_r/mask.mrc:160:open ex_cdisc_A:$PB/cdisc_A:$CD/node_r/mask.mrc:225.456:A \
  ex_cdisc_closed:$PB/cdisc_closed:$CD/node_r/mask.mrc:225.456:closed > $OUT/ex_switches/maptests.txt 2>&1 ) &
newdir cahra_open_bc
{ $PY $T/pick_score_bc.py score.txt open cahra $GT $CA > ../pick.txt 2>&1 && $PY $T/motion_bins4.py score.txt states.txt cahra $GT tree $CA > ../bins.txt 2>&1 \
  && recon cahra_open_bc $CA/node_r/stageB/b.simple $CA/node_r/recB0/hetero_step.sh $CA/node_r/mask.mrc open; } || log "cahra_open_bc FAILED: $(tail -qn1 ../pick.txt ../bins.txt 2>&1)"
for ex in tree truth; do newdir cdisc_A_${ex}junk; cp $PB/cdisc_A/work/score.txt .
  { $PY $T/motion_bins4.py score.txt states.txt cdisc $ST $ex $CD > ../bins.txt 2>&1 \
    && recon cdisc_A_${ex}junk $CD/node_r/stageB/b.simple $CD/node_r/recB0/hetero_step.sh $CD/node_r/mask.mrc A; } || log "cdisc_A_${ex}junk FAILED: $(tail -1 ../bins.txt)"
done
finish

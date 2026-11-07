#!/bin/bash
#SBATCH -J klem_calE_igm
#SBATCH -p csbdevel
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 00:28:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Swing calibration of the switch / swing call: for every IgM arm angle the M tree's best population-wide kept-site posterior
# (pick_igm.py), 5 at-node bins along it (motion_bins4.py: tree junk out, no truth in the edges), the root's recB0 reconstruct3D
# with nstates=5, postprocess=no (seriesE_lib.sh); phaseE_maptests.py on each and on the three phase-C IgM series
set -uo pipefail
NP=2; NT=6; TAG=calE_igm; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/igm; mkdir -p $OUT/scores; IG=$S/fulltoy_m52M_igm/1_hetero_tree; IT=$K/factorial/armsigm1_truth/arms_truth.txt
I=$K/factorial/motion_series_igm; MI=$IG/node_r/mask.mrc; MB=$S/fulltoy_igm/1_hetero_tree/node_r/mask.mrc
for e in m52M_r_s1_a4th_q5:$MI m52M_r2_s4_a2th_q5:$MI m8c17_r2_s1_a2th_q5:$MB; do IFS=: read nm mk <<< "$e"
  mkdir -p $OUT/ex_$nm; ( $PY $T/phaseE_maptests.py ex_$nm:$I/$nm:$mk:220:mean > $OUT/ex_$nm/maptests.txt 2>&1 ) &
done
$PY $T/pick_igm.py $OUT/scores $IG $IT > $OUT/picks.txt 2>&1 || { log "pick FAILED: $(tail -2 $OUT/picks.txt)"; finish; exit 1; }
mapfile -t PK < <(grep "^PICK" $OUT/picks.txt | grep -v " same$"); log "${#PK[@]} distinct picks"
for l in "${PK[@]}"; do read _ ang site j rho f <<< "$l"; nm=m52M_${site}_p${j}_${ang/_/}; newdir $nm
  $PY $T/motion_bins4.py $f states.txt arms $IT tree $IG $ang > ../bins.txt 2>&1 || { log "$nm bins FAILED: $(tail -1 ../bins.txt)"; continue; }
  recon $nm $IG/node_r/stageB/b.simple $IG/node_r/recB0/hetero_step.sh $MI mean
done
finish

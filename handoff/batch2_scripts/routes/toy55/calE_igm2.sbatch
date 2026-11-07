#!/bin/bash
#SBATCH -J klem_calE_igm2
#SBATCH -p csbdevel
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 00:50:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Swing calibration on a second IgM tree, calE_igm.sbatch's recipe on the m8c17 baseline (stree_port/fulltoy_igm): every arm angle's
# best population-wide kept-site posterior (pick_igm.py), 5 at-node bins (tree junk out), the root's recB0 reconstruct3D with
# nstates=5, postprocess=no; phaseE_maptests.py on each
set -uo pipefail
NP=2; NT=6; TAG=calE_igm2; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/igm_m8c17; mkdir -p $OUT/scores; IG=$S/fulltoy_igm/1_hetero_tree; IT=$K/factorial/armsigm1_truth/arms_truth.txt
$PY $T/pick_igm.py $OUT/scores $IG $IT > $OUT/picks.txt 2>&1 || { log "pick FAILED: $(tail -2 $OUT/picks.txt)"; exit 1; }
mapfile -t PK < <(grep "^PICK" $OUT/picks.txt | grep -v " same$"); log "${#PK[@]} distinct picks"
for l in "${PK[@]}"; do read _ ang site j rho f <<< "$l"; nm=m8c17_${site}_p${j}_${ang/_/}; newdir $nm
  $PY $T/motion_bins4.py $f states.txt arms $IT tree $IG $ang > ../bins.txt 2>&1 || { log "$nm bins FAILED: $(tail -1 ../bins.txt)"; continue; }
  recon $nm $IG/node_r/stageB/b.simple $IG/node_r/recB0/hetero_step.sh $IG/node_r/mask.mrc mean
done
finish

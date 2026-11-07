#!/bin/bash
#SBATCH -J klem_calE_igm3
#SBATCH -p csbdevel
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 00:35:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Swing calibration at the sites' own nodes on the m8c17 IgM baseline (stree_port/fulltoy_igm), whose kept sites are scored on their
# own levels: every arm angle's best kept-site posterior on a level of >= 20k rows (pick_igm2.py), 5 at-node bins over the level
# (tree junk out), the level node's own project, recB0 reconstruct3D with nstates=5, postprocess=no, and mask; phaseE on each
set -uo pipefail
NP=2; NT=6; TAG=calE_igm3; source /mnt/beegfs/elmlund/afan/routes/toy55/seriesE_lib.sh
OUT=$K/factorial/phaseE/igm_m8c17; mkdir -p $OUT/scores; IG=$S/fulltoy_igm/1_hetero_tree; IT=$K/factorial/armsigm1_truth/arms_truth.txt
$PY $T/pick_igm2.py $OUT/scores $IG $IT 20000 > $OUT/picks2.txt 2>&1 || { log "pick FAILED: $(tail -2 $OUT/picks2.txt)"; exit 1; }
mapfile -t PK < <(grep "^PICK" $OUT/picks2.txt | grep -v " same$"); log "${#PK[@]} distinct picks"
for l in "${PK[@]}"; do read _ ang site j rho f lev nr <<< "$l"; nm=m8c17_${site}_p${j}_${ang/_/}; N=$IG/node_$lev; newdir $nm
  [ -s $N/n.simple ] && [ -s $N/recB0/hetero_step.sh ] || { log "$nm: node_$lev lacks n.simple or recB0"; continue; }
  $PY $T/motion_bins4.py $f states.txt arms $IT tree $IG $ang > ../bins.txt 2>&1 || { log "$nm bins FAILED: $(tail -1 ../bins.txt)"; continue; }
  recon $nm $N/n.simple $N/recB0/hetero_step.sh $N/mask.mrc mean
done
finish

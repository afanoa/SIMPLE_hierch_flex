#!/bin/bash
#SBATCH -J klem_groupdiscrim
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# what separates the merged groups of the root split test (env MODE, T (splitparent_k6), V (NUIS 0|1), G1 G2 = the groups to compare)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; T=${T:-splitparent_k6}; V=${V:-0}
RA=$A/rnasep_klem/synth_${MODE}_A
for H in A B; do R=$A/rnasep_klem/synth_${MODE}_$H; D=$R/fam0/region$T/level1
  $PY $W/group_labels.py $D/states.txt $RA/fam0/region$T/merge_nuis${V}_$H.npz $D/groups_nuis$V.txt 2>&1 | grep -v Warning
  echo "== $MODE half $H: group ${G1:-1} vs group ${G2:-2} (NUIS=$V)"
  $PY $W/split_discrim.py $W/sim${MODE}_$H/all.txt $W/lat_${MODE}_$H/ptcl3D_r0.txt $D/groups_nuis$V.txt:${G1:-1} $D/groups_nuis$V.txt:${G2:-2} $W/sim${MODE}_$H/gt_particles.txt 2>&1 | grep -v Warning | head -14; done

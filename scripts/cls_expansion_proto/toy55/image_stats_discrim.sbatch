#!/bin/bash
#SBATCH -J klem_imgstats
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# image statistics of the merged root groups (env MODE, T split dir tag, V nuis variant + band, G1 G2)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; T=${T:-splitparent_k6}; V=${V:-nuis0}
RA=$A/rnasep_klem/synth_${MODE}_A
for H in A B; do R=$A/rnasep_klem/synth_${MODE}_$H; D=$R/fam0/region$T/level1; NPZ=$(ls $RA/fam0/region$T/merge_${V}_$H.npz $RA/fam0/region$T/merge_${V}_b40_$H.npz 2>/dev/null | head -1)
  $PY $W/group_labels.py $D/states.txt $NPZ $D/groups_$V.txt 2>&1 | grep -v Warning
  echo "== $MODE half $H ($T, $V): image statistics, group ${G1:-1} vs ${G2:-2}"; $PY $W/image_stats_discrim.py $W/sim${MODE}_$H/all.mrcs $D/groups_$V.txt ${G1:-1} ${G2:-2} $W/sim${MODE}_$H/gt_particles.txt 2>&1 | grep -v Warning; done

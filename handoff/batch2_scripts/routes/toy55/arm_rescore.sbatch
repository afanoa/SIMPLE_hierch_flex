#!/bin/bash
#SBATCH -J klem_armrescore
#SBATCH -p norm
#SBATCH --qos=norm
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# a truth coordinate every particle has: the closed-arm occupancy of the map it was simulated from (arm_readout.py on every
# conformer map of the set), then the fast run's factors re-scored with it. env TAGS (space list), EXTRA (codes per tag via env CODES_<tag>)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; export OMP_NUM_THREADS=8
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for tag in ${TAGS:?}; do
  tab=$K/armtab_$tag.txt
  if [ ! -s $tab ]; then
    $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $K/armtab_$tag.raw $T/sim$tag/maps/conf*_m160.mrc > /dev/null 2>&1
    awk '!/^#/ {if (match($1, /conf[0-9]+_/)) print substr($1, RSTART + 4, RLENGTH - 5), $2}' $K/armtab_$tag.raw > $tab
  fi
  echo "== $tag: closed-arm occupancy of $(grep -c . $tab) conformer maps (range $(sort -k2 -g $tab | head -1 | cut -d' ' -f2) .. $(sort -k2 -g $tab | tail -1 | cut -d' ' -f2))"
  case $tag in swingbreath*) codes=swing:200:1.0,breath:300:1.0;; prod*) codes=grid:400:2.0:2.0;; *) codes=swing:200:1.0;; esac
  nf=$(cat $K/ff2_$tag/expand/nfam.txt)
  $PY $T/factor_sample.py $K/ff2_$tag $K/ff2_$tag/factors_arm --fams $(seq -s, 1 $nf) --glob 'region_c*' --nbin 8 --budget 1000 --minend 300 --gt $T/sim$tag/gt_particles.txt --codes $codes,table:$tab:closedarm 2>&1 | grep -v Warning | grep 'factors at\|factor [0-9]\|closedarm\|truth swing\|truth grid_k'
done

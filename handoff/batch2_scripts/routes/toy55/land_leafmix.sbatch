#!/bin/bash
#SBATCH -J land_leafmix
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_leafmix.%j.out
# whole-set map-level populations: the 64 plain f692st leaf half-maps unmixed onto the deposited conformers, pooled by particle counts
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; U=$K_/tw_units/f692st_plain; O=$K_/landscape/leafmix_f692st; mkdir -p $O
CNT=$(grep -v "^#" $U/units_list.txt | awk '{print $4}' | paste -sd, -)
for hv in even odd; do
  $PY $T/land_confmix.py $T/maps160/toy_consensus.mrc $K_/landscape/confs/arm_confs.txt $O/confmix_$hv.txt --refs $T/maps160/docked_all/conf*_m160.mrc $T/maps160/holo_docked/conf*_m160.mrc \
     --maps $(grep -v "^#" $U/units_list.txt | awk -v c=$([ $hv = even ] && echo 2 || echo 3) '{print $c}') --counts $CNT 2>&1 | grep --line-buffered -v Warning | sed "s/^/$hv /"
done
echo "[$(date '+%F %T')] land_leafmix DONE"

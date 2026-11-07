#!/bin/bash
#SBATCH -J null_xem
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# no-heterogeneity null of the cross-half final pass: one conformer's particles in NK random groups per half -> the same
# split test + cross-half EM + re-test as tree_finalize (env MODE, NK 4, CONF, POSES consensus|truth); K* > 1 = false states
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; MODE=${MODE:-55_700k}
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16; O=$K/null_xem_${MODE}_${NK:-4}${TAGX:-}; mkdir -p $O
$PY $W/null_labels.py $W/sim${MODE}_A/gt_particles.txt $W/sim${MODE}_B/gt_particles.txt $O/cand_A.txt $O/cand_B.txt ${NK:-4} ${CONF:-} || exit 1
PA=$W/lat_${MODE}_A/ptcl3D_r0.txt; PB=$W/lat_${MODE}_B/ptcl3D_r0.txt; [ "${POSES:-consensus}" = truth ] && { PA=$W/sim${MODE}_A/all.txt; PB=$W/sim${MODE}_B/all.txt; }
$PY -u $W/splitbench7.py --sim $W/sim${MODE}_A,$W/sim${MODE}_B --poses $PA,$PB --cand $O/cand_A.txt,$O/cand_B.txt \
  --gt $W/sim${MODE}_A/gt_particles.txt,$W/sim${MODE}_B/gt_particles.txt --out $O --cache $W/benchcache_$MODE --xem ${XEM:-6} --em_lim ${EMLIM:-20} --em_min ${EMMIN:-1000} --force_em 1 --retest_beyond_em 1 --pair hungarian ${EXTRA:-} 2>&1 | grep --line-buffered -v Warning

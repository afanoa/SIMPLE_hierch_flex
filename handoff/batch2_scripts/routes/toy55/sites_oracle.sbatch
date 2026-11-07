#!/bin/bash
#SBATCH -J klem_sites_oracle
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# simsites per-site truth table (site_truth.py) + the ORACLE per-particle ceiling of every site (oracle_sites.py: true combination maps,
# per-image scale, whitened, independent coefficients) at the TRUE poses and at the REFINED poses, lp 20 / 15 / 10. Outputs rnasep_klem/factorial/simsites_truth
set -uo pipefail
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=$K_/factorial/simsites_truth; mkdir -p $O
[ -s $O/sites_truth.txt ] || $PY $T/site_truth.py sites $T/simsites/gt_particles.txt $T/maps160/simsites/combos.txt $O/sites_truth.txt || exit 1
for pz in ${POSESETS:-true refined}; do P=$T/simsites/all.txt; [ $pz = refined ] && P=$T/lat_sites/ptcl3D_r0.txt
  echo "== oracle at the $pz poses ($P)"; $PY -u $T/oracle_sites.py --sim $T/simsites --poses $P --truth $O/sites_truth.txt --maps $T/maps160/simsites --out $O/oracle_$pz.npz --nsub ${NSUB:-40000} --lp 20,15,10 2>&1 | grep -v Warning; done

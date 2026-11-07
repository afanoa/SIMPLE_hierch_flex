#!/bin/bash
#SBATCH -J toy10_levers_density
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Lever 2: states as modes of the deconvolved latent density (density_modes.py) on existing toy boards: the parent
# (klem_junk), the feedback child boards (klem_fb_fam1/2) and any boards named in DIRS; readouts vs truth against
# k-means K and the argmax cells. env: DIRS (space-separated klem dirs under lat/), K0 (20), KM (10)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W
for d in ${DIRS:-klem_junk klem_fb_fam1 klem_fb_fam2}; do
  echo "== $d"; $PY $W/density_modes.py $W/lat/$d --gt $W/sim/gt_particles.txt --k0 ${K0:-20} --kmeans ${KM:-10} 2>&1 | grep -v Warning
done
echo "[$(date '+%F %T')] density DONE"

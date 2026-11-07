#!/bin/bash
#SBATCH -J toy10_klem_zoom
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Hierarchical zoom on a toy: NB bins along the main axis of a finished klem latent (KDIR), then a fresh klem run on the
# particles of each bin alone (same 7x7 grid, PCA init): does purity within a bin exceed the global run's? env: TAG, KDIR, NB (3), INIT (pca)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; TAG=${TAG:-}; OUT=$W/lat$TAG; KDIR=${KDIR:-$OUT/klem_junk}; NB=${NB:-3}; Z=$OUT/zoom; mkdir -p $Z; cd $W
NPTCL=$(grep -c . $W/${SIMDIR:-sim}/gt_particles.txt); $PY $W/klem_bins.py $KDIR none $Z/states.txt --nb $NB --nptcl $NPTCL
for b in $(seq 1 $NB); do awk -v b=$b '$1 == b {print NR}' $Z/states.txt > $Z/subset_bin$b.txt; echo "bin $b: $(wc -l < $Z/subset_bin$b.txt) particles"; done
for b in $(seq 1 $NB); do echo "[$(date '+%F %T')] zoom bin $b"; env TAG=$TAG KTAG=_zoom$b INIT=${INIT:-pca} NSUB=0 SUBSET=$Z/subset_bin$b.txt bash $W/run_klem.sbatch 2>&1 | grep "final:\|subset\|particles (\|per-conformer mean" | sed "s/^/[bin $b] /"; done
echo "[$(date '+%F %T')] zoom DONE"

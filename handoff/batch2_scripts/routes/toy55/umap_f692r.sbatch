#!/bin/bash
#SBATCH -J klem_umap
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=200G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-particle UMAP of the real RNase P states of a consolidation (env C = dir with states.txt over master rows, arm_states.txt; TAG):
# master-row labels -> per-half labels (meta_rows), one cross-half E-step over the states at the family-refined poses (state_graph.py,
# |CTF|), UMAP of the memberships (umap_states.py) coloured by state and by the closed-arm readout of each particle's state.
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; C=${C:?}; TAG=${TAG:?}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8; P=$K/graph_$TAG; S=$W/simrnasep692; PO=$K/dtree_simrnasep692_f692r/poses_famrefined
for h in A B; do $PY -c "import numpy as np; st = np.loadtxt('$C/states.txt', dtype=int); mr = np.loadtxt('${S}_$h/meta_rows.txt').astype(int) - 1; np.savetxt('${P}_labels_$h.txt', st[mr], fmt='%d')"; done
echo "labels: $(awk '$1 > 0' ${P}_labels_A.txt | wc -l) / $(awk '$1 > 0' ${P}_labels_B.txt | wc -l) labelled particles per half"
$PY $W/state_graph.py --sim ${S}_A,${S}_B --poses ${PO}_A.txt,${PO}_B.txt --labels ${P}_labels_A.txt,${P}_labels_B.txt --cache $W/benchcache_rnasep692 \
  --png $P.png --ctfabs 1 --title "f692r $TAG" 2>&1 | grep --line-buffered -v Warning | tail -12
$PY $W/umap_states.py $P.npz ${P}_umap.png --arm $C/arm_states.txt --nsub ${NSUB:-80000} --title "f692r $TAG: per-particle UMAP of the cross-half state memberships" 2>&1 | grep -v Warning | tail -3

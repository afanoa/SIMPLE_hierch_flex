#!/bin/bash
#SBATCH -J klem_coord
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=120G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# stage 4 of the general pipeline: the coordinate method (coord_loop.sbatch) inside the merge-up's final states and branches, where
# hard labels cannot split a continuum. A branch node (>= 3 final states) starts from its states; a single final state starts from a
# random SEEDK split of its particles. Nodes run one after another in this job. Output per node: OUT/<name>/final_coord.txt (the
# per-particle coordinate), final_states.txt (quantile bins, maps in final/), final_even_bins.txt (the even-width sampled trajectory).
# env TAG, UNITS (per-row final units), BASE (project with the refined poses), NODES ('name:u+u+u;name:u' or auto = every final state
# with >= MINP particles), MINP (15000), SEEDK (8), KL (8), NIT (2), OUT, CODES (truth codes, environment)
set -uo pipefail
T=/mnt/beegfs/elmlund/afan/routes/toy55; PY=/mnt/nasapps/production/python/3.12/bin/python3; O=${OUT:?}; mkdir -p $O; cd $O
log(){ echo "[$(date '+%F %T') coord-nodes ${TAG:?}] $*"; }; U=${UNITS:?}; NL=${NODES:-auto}
if [ "$NL" = auto ]; then NL=$(awk -v m=${MINP:-15000} '$1 > 0 {c[$1]++} END {for (u in c) if (c[u] >= m) printf "state%s:%s;", u, u}' $U); fi
log "nodes: $NL"
for spec in $(echo "$NL" | tr ';' ' '); do nm=${spec%%:*}; us=${spec#*:}; [ -s $O/$nm/final_coord.txt ] && { log "$nm done before"; continue; }
  $PY - <<PYEOF
import numpy as np
u = np.loadtxt('$U', dtype=int); us = [int(x) for x in '$us'.split('+')]; lab = np.zeros(len(u), int); rng = np.random.default_rng(7)
if len(us) >= 3:
    for i, x in enumerate(us, 1): lab[u == x] = i
else:
    m = np.isin(u, us); lab[m] = rng.integers(1, ${SEEDK:-8} + 1, m.sum())
np.savetxt('$O/${nm}_lab0.txt', lab, fmt='%d'); print('node $nm: units $us, %d particles, %d starting labels' % ((lab > 0).sum(), len(np.unique(lab[lab > 0]))))
PYEOF
  log "coordinate loop on node $nm"; TAG=$TAG OUT=$O/$nm BASE=${BASE:?} LAB0=$O/${nm}_lab0.txt KL=${KL:-8} NIT=${NIT:-2} bash $T/coord_loop.sbatch 2>&1 | grep -a "TRUTH\|PCA variance\|even sampling\|DONE\|FAILED" | sed "s/^/   $nm: /" | cut -c1-260
done; log DONE

#!/bin/bash
#SBATCH -J klem_lsamp_chain
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=220G
#SBATCH -t 16:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape sampling of a finished tree in one job: (1) score space with the pixel split + per-particle noise levels (lsamp_score2.sbatch's
# steps), (2) every terminal node (lsamp_tree.sbatch's steps: lsamp2.py with the defaults + reconstruct3D + readouts).
# env TREE (dir with leaves_list.txt + final_labels.txt), OUT (score dir; nodes in OUT/tree), SIM, MR, CP, CTFABS, MSKDIAM, POSES (optional),
# ARGS (extra lsamp2.py flags), NODES (optional), VIEWBASE, GT (synthetic)
set -uo pipefail
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; log(){ echo "[$(date '+%F %T') lsamp_chain $(basename ${OUT:?})] $*"; }
log "scoring $TREE"
PY=/mnt/nasapps/production/python/3.12/bin/python3
if [ ! -s $OUT/meta.npz ] || ! $PY -c "import numpy as np, sys; sys.exit(0 if 'ppnoise' in np.load('$OUT/meta.npz').files else 1)" 2>/dev/null; then
  TREE=$TREE OUT=$OUT SIM=$SIM MR=$MR CP=$CP CTFABS=${CTFABS:-0} EMLIM=${EMLIM:-12} MSKDIAM=${MSKDIAM:-230} POSES=${POSES:-} bash $T/lsamp_score2.sbatch; fi
[ -s $OUT/meta.npz ] || { log "scoring FAILED"; exit 1; }
P=${POSES:-$OUT/poses_combined.txt}
log "terminal nodes"
TREE=$TREE SC=$OUT OUT=$OUT/tree SIM=$SIM MR=$MR CP=$CP CTFABS=${CTFABS:-0} MSKDIAM=${MSKDIAM:-230} POSES=$P ARGS="${ARGS:---ends centre}" NODES="${NODES:-}" \
  VIEWBASE=${VIEWBASE:-$A/rnasep_klem/view_rflex/lsamp_$(basename $OUT)} GT=${GT:-} bash $T/lsamp_tree.sbatch
log "DONE"

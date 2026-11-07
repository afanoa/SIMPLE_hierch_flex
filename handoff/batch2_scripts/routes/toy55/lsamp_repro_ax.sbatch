#!/bin/bash
#SBATCH -J klem_lsamp_reprox
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# lsamp_repro.sbatch with the stage options in LSARGS (default --axes 2: the second score axis as its own 1-D landscape): (1) score space against the run's leaves (lsamp_score2.sbatch's steps,
# pixel split + per-particle noise levels) unless SC has it, (2) lsamp2h.py with the defaults (auto axes + check-3 gate, centred ends,
# structural spacing) -> samples, coords.npz, held-out / own-pixel / final crop maps, (3) arm readout of every map set.
# env RUN, TREE (leaves_list.txt + final_labels.txt), SC (score dir), PREFIX (leaf-name prefix of the node, e.g. r1_), OUT, POSES (optional), LSARGS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${OUT:?}; mkdir -p $O
log(){ echo "[$(date '+%F %T') lsamp_repro ${RUN:?}] $*"; }; H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692
if [ ! -s ${SC:?}/meta.npz ] || ! $PY -c "import numpy as np, sys; sys.exit(0 if 'ppnoise' in np.load('$SC/meta.npz').files else 1)" 2>/dev/null; then
  log "scoring $TREE"; TREE=${TREE:?} OUT=$SC SIM=${H}_A,${H}_B MR=$MR CP=$CP CTFABS=1 EMLIM=12 MSKDIAM=232 POSES=${POSES:-} bash $T/lsamp_score2.sbatch 2>&1 | grep -v "dumped"; fi
P=${POSES:-$SC/poses_combined.txt}; log "landscape of the ${PREFIX:?} leaves"
$PY -u $T/lsamp2h.py --scores $SC --sim ${H}_A,${H}_B --mr $MR --cache $CP --poses $P --ctfabs 1 --mskdiam 232 --out $O --use $PREFIX ${LSARGS:---axes 2} --ends centre --tag $RUN 2>&1 \
  | grep --line-buffered -v "Warning\|warn(" > $O/lsamp.log
grep "sampled axes\|dropping\|resolvable samples\|CHOSEN\|steps CV\|maps (cells" $O/lsamp.log | cut -c1-200
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
for d in heldout_A heldout_B final_all; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm_$d.txt $O/$d 2>&1 | grep -v Warning > /dev/null
  echo "$d: $(awk '!/^#/ {k = $1; sub(/.*state_/, "", k); sub(/\.mrc/, "", k); printf "%s %.2f/%.2f  ", k + 0, $2, $3}' $O/arm_$d.txt)"; done
log DONE

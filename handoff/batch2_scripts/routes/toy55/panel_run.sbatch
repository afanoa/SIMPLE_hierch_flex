#!/bin/bash
#SBATCH -J klem_panel
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# one benchmark set through the general pipeline, frozen settings: (1) discovery on all particles = K-state refine3D_states from random
# seeds (pose_policy local, 20 -> 12 A, 6 its); (2) leaf tree from the leaves' half-maps (state_merge.py, complete linkage, 40-10 A);
# (3) merge-up along the tree (bu_merge.sbatch: splitbench8 on the even/odd halves, label-shuffle null); scorecard after every stage.
# env TAG (sim<TAG>, proj<TAG>, lat_<TAG>), K (20), NSUB (all rows), OUT (rnasep_klem/panel_<TAG>); scorecard GT/CLASSES/CODES (environment)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; TAG=${TAG:?}; K=${K:-20}
O=${OUT:-$K_/panel_$TAG}; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') panel $TAG] $*"; }; GT=${GT:-$W/sim$TAG/gt_particles.txt}
NR=$(grep -c . $W/lat_$TAG/ptcl3D_r0.txt)
[ -s seed_states.txt ] || $PY -c "import numpy as np; rng = np.random.default_rng(7); n = $NR; s = np.zeros(n, int); r = rng.choice(n, min(n, ${NSUB:-$NR}), replace=False); s[r] = rng.integers(1, $K + 1, len(r)); np.savetxt('seed_states.txt', s, fmt='%d')"
log "discovery: $K states from random seeds over $(awk '$1 > 0' seed_states.txt | wc -l) of $NR particles"
if [ ! -s refined_states.txt ]; then
  cp -f $W/proj$TAG/proj$TAG.simple p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
  simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local pgrp=c1 mskdiam=220 lpstart=20 lpstop=12 maxits=6 ml_reg=no nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR' ref.log | cut -c1-150)"; exit 1; }
  RD=$(ls -td $O/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > refined_states.txt
fi
RD=$(ls -td $O/*_refine3D_states | head -1); sc(){ $PY $W/scorecard.py $1 $GT ${CLASSES:+--classes $CLASSES} ${CODES:+--codes $CODES} --tag "$2" 2>&1 | grep -v Warning; }
sc refined_states.txt "discovery leaves"
[ -s smerge/merge_tree.txt ] || $PY $W/state_merge.py smerge --rdir $RD --labels refined_states.txt --gt $GT --unfil --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "^== \|leaf pair" | cut -c1-200
log "merge-up along the leaf tree"; TAG=$TAG LEAVES=$O/refined_states.txt TREE=$O/smerge/merge_tree.txt PROJ=$RD/p.simple OUT=$O/bu GT=$GT PAR=${PAR:-8} bash $W/bu_merge.sbatch 2>&1 | grep -v "^\[.*eo:" | tail -40
log DONE

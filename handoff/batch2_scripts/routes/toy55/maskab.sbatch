#!/bin/bash
#SBATCH -J maskab
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=96G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Focus-mask A/B at one rflex_s30a node: rflex3's focused K-state refinement repeated EXACTLY (node project n.simple, selection with
# the node's flex_states.txt, the same flex-state refs, refine3D prob_neigh geom lp 10 maxits 8 pcg_focus nspace 5000 box_crop 88,
# nparts 4 x nthr 16, focus3 snapshot) with ONLY pcg_mskfile = MASK (default: the cheap children's focus mask). Then truth
# (scorecard.py) and per-state map scores (subcls_statescore.py) for this run and the node's current-mask refinement, and the
# agreement of the two labellings. Nothing in the rflex_s30a node dir is written. env NODE (r12 | r22), MASK, TAG (swing30)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
NODE=${NODE:?}; TAG=${TAG:-swing30}; GT=$T/sim$TAG/gt_particles.txt; TM=$T/sim$TAG/maps; MD=220
D=$K_/rflex_s30a/node_$NODE; MASK=${MASK:-$K_/cheapkids_s30a_$NODE/kids_mask.mrc}; O=$K_/maskab_s30a_${NODE}_$(basename $MASK .mrc); mkdir -p $O && cd $O
log(){ echo "[$(date '+%F %T') maskab $NODE] $*"; }
[ -s $MASK ] && [ -s $D/n.simple ] && [ -s $D/flex_states.txt ] || { log "inputs missing"; exit 1; }
KS=$(grep -c . $D/members.txt); VOLS=""; for j in $(seq 1 $KS); do VOLS="$VOLS vol$j=$D/refs/ref_state$(printf %02d $j).mrc"; done
log "mask $MASK; current mask $D/mask.mrc; $KS states from $D/refs"
t0=$(date +%s)
if [ -z "$(ls $O/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cp -f $D/n.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$D/flex_states.txt mkdir=no > rsel.log 2>&1 || { log "state selection FAILED"; exit 1; }
  simple_exec prg=refine3D projfile=r.simple $VOLS nstates=$KS refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=10 maxits=8 center=no rec_backend=pcg pcg_mskfile=$MASK pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }
fi
RD=$(ls -td $O/*_refine3D | head -1); RP=$RD/r.simple; [ -s $RP ] || RP=$O/r.simple; t1=$(date +%s); log "refine3D $((t1 - t0)) s, $(grep -a -c '>>> ITERATION' ref.log) iterations"
simple_exec prg=print_project_field projfile=$RP oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/refined_states.txt
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$T'); a = np.loadtxt('$D/refined_states.txt').astype(int); b = np.loadtxt('$O/refined_states.txt').astype(int); m = (a > 0) & (b > 0)
_, ai = np.unique(a[m], return_inverse=True); _, bi = np.unique(b[m], return_inverse=True); J = np.zeros((ai.max() + 1, bi.max() + 1)); np.add.at(J, (ai, bi), 1); J /= J.sum()
pa, pb = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(pa, pb)[nz])).sum(); h = lambda p: -(p[p > 0] * np.log(p[p > 0])).sum()
print('labelling agreement current-mask vs new-mask refinement: NMI %.3f; same state for %.0f %% of particles (state numbers follow the shared refs)' % (I / np.sqrt(h(pa) * h(pb)), 100 * (a[m] == b[m]).mean()))
PYEOF
echo "== truth scorecards (node $NODE members)"
$PY $T/scorecard.py $D/flex_states.txt $GT --codes swing:200:1 --tag "flex states (start)" 2>&1 | grep -v Warning
$PY $T/scorecard.py $D/refined_states.txt $GT --codes swing:200:1 --tag "refined, CURRENT mask" 2>&1 | grep -v Warning
$PY $T/scorecard.py $O/refined_states.txt $GT --codes swing:200:1 --tag "refined, $(basename $MASK .mrc)" 2>&1 | grep -v Warning
echo "== per-state maps vs truth"
$PY $T/subcls_statescore.py $(ls -td $D/*_refine3D | head -1) $D/refined_states.txt $GT $TM --mskdiam $MD --tag "CURRENT mask" 2>&1 | grep -v Warning
$PY $T/subcls_statescore.py $RD $O/refined_states.txt $GT $TM --mskdiam $MD --tag "$(basename $MASK .mrc)" 2>&1 | grep -v Warning
log "DONE ($(( $(date +%s) - t0 )) s)"

#!/bin/bash
#SBATCH -J klem_lsamp_disc
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=180G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# leaves for validating landscape sampling on a synthetic set without a tree: K-state refine3D_states from random seeds (panel_run.sbatch's
# discovery step: pose_policy local, 20 -> 12 A, 6 its) + the complete-linkage tree over the states (state_merge.py). Writes the rflex-like
# layout lsamp_score.sbatch reads: OUT/leaves_list.txt (name even odd count), OUT/final_labels.txt, OUT/final_tree/merge_tree.txt.
# env TAG (proj<TAG>, sim<TAG>), K (16), OUT (rnasep_klem/lsamp/disc_<TAG>), MSKDIAM (220)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=12
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; TAG=${TAG:?}; K=${K:-16}
O=${OUT:-$K_/lsamp/disc_$TAG}; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') lsamp_disc $TAG] $*"; }; GT=$T/sim$TAG/gt_particles.txt
NR=$(simple_exec prg=print_project_field projfile=$T/proj$TAG/proj$TAG.simple oritype=ptcl3D mkdir=no 2>/dev/null | grep -c .)
[ -s seed_states.txt ] || $PY -c "import numpy as np; rng = np.random.default_rng(7); np.savetxt('seed_states.txt', rng.integers(1, $K + 1, $NR), fmt='%d')"
log "discovery: $K states from random seeds over $NR particles"
if [ ! -s refined_states.txt ]; then
  cp -f $T/proj$TAG/proj$TAG.simple p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=12 mkdir=no > up.log 2>&1
  simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local pgrp=c1 mskdiam=${MSKDIAM:-220} lpstart=20 lpstop=12 maxits=6 ml_reg=no nparts=4 nthr=12 > ref.log 2>&1 \
    || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR' ref.log | cut -c1-150)"; exit 1; }
  RD=$(ls -td $O/*_refine3D_states | head -1)
  simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > refined_states.txt
fi
RD=$(ls -td $O/*_refine3D_states | head -1); : > leaves_list.txt
for s in $(seq 1 $K); do n=$(awk -v s=$s '$1 == s' refined_states.txt | wc -l); b=$RD/recvol_state$(printf %02d $s); [ $n -gt 0 ] && echo "s$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc $n" >> leaves_list.txt; done
awk 'NR == FNR {sub(/^s0*/, "", $1); id[$1] = FNR; next} {print ($1 in id) ? id[$1] : 0}' leaves_list.txt refined_states.txt > final_labels.txt   # leaf index (leaves_list order)
log "leaves: $(wc -l < leaves_list.txt) states, sizes $(awk '{printf "%s ", $4}' leaves_list.txt)"
[ -s final_tree/merge_tree.txt ] || $PY $T/state_merge.py final_tree --rdir $RD --labels refined_states.txt --gt $GT --unfil --link complete --ncut 8,4,2 2>&1 | grep -v Warning | grep "^== \|leaf pair" | cut -c1-200
$PY $T/scorecard.py refined_states.txt $GT --tag "discovery states" 2>&1 | grep -v Warning | head -5
log DONE

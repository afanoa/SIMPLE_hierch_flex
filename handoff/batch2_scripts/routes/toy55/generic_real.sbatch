#!/bin/bash
#SBATCH -J klem_generic
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=240G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the general pipeline end to end on real RNase P (all particles, no family split), frozen settings:
# (0) seeds = the whole-particle board's cells in K spectral clusters (BOARD: an existing generic board, e.g. f692r/parent);
# (1) discovery: K-state refine3D_states on all particles from those seeds (local poses, 20 -> 8 A, 6 its, distributed NP x 8);
# (2) leaf tree from the leaves' half-maps (state_merge.py, complete linkage); (3) merge-up along the tree on the real half-stacks
# (bu_merge_real.sbatch: splitbench9, cross-half EM, stratified null, |CTF|); (4) final states -> one cross-half scoring pass ->
# per-node axes, multi-axis even sampling (coord_tree2.py) -> (5) the sampled cells reconstructed + closed-arm readout (recon_nodes).
# env RTAG (f692g), K (32), SRC (master project), BOARD (f692r/parent), NP (20), MSKDIAM (232); SNAP (focus3: > vol20)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RTAG:-f692g}; R=$K_/$RT
K=${K:-32}; NP=${NP:-20}; MD=${MSKDIAM:-232}; SRC=${SRC:-$A/rnasep_coupled/full692k_2d/rnasep.simple}; BD=${BOARD:-$K_/f692r/parent}
mkdir -p $R/disc && cd $R; log(){ echo "[$(date '+%F %T') generic $RT] $*"; }; NR=$(simple_exec prg=print_project_field projfile=$SRC oritype=ptcl3D mkdir=no 2>/dev/null | grep -c .)
log "master project $SRC: $NR rows; board $BD; K $K"
# ---- (0) seeds
[ -s $R/seed_states.txt ] || $PY - <<PYEOF 2>&1 | grep -v Warning
import glob, numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23])).reshape(n, n, n).astype(float)
fs = sorted(glob.glob('$BD/cells/cell_*.mrc')); V = np.array([rd(f) for f in fs]); D = V - V.mean(0); X = np.array([np.fft.fftn(d).ravel() for d in D]); X /= np.linalg.norm(X, axis=1)[:, None]
S = np.maximum(np.real(X @ X.conj().T), 0); np.fill_diagonal(S, 0); dg = S.sum(1) + 1e-9; ev, U = np.linalg.eigh(np.eye(len(S)) - S / np.sqrt(np.outer(dg, dg))); k = min($K, len(fs))
E = U[:, :k] / np.maximum(np.linalg.norm(U[:, :k], axis=1)[:, None], 1e-12); rng = np.random.default_rng(0); best = None
for rep in range(20):
    C = E[rng.choice(len(E), k, replace=False)]
    for it in range(100): lab = np.argmin(((E[:, None] - C[None]) ** 2).sum(-1), 1); C = np.array([E[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
    sse = ((E - C[lab]) ** 2).sum(); best = (sse, lab) if best is None or sse < best[0] else best
lab = best[1]; L_ = np.loadtxt('$BD/klem_latent.txt', comments='#'); pind, cell = L_[:, 0].astype(int), L_[:, 3].astype(int); ncell = len(fs)
out = np.zeros($NR, int); on = (cell >= 1) & (cell <= ncell); out[pind[on] - 1] = lab[cell[on] - 1] + 1
np.savetxt('$R/seed_states.txt', out, fmt='%d'); print('board seeds: %d cells in %d spectral clusters, %d particles seeded' % (ncell, k, (out > 0).sum()))
PYEOF
# ---- (1) discovery: K-state refine3D_states on all particles (distributed)
cd $R/disc; if [ ! -s $R/leaves.txt ]; then cp -f $SRC p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$R/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
  export SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8; KS=$(awk '$1 > 0' $R/seed_states.txt | sort -n | uniq | wc -l)
  log "discovery: refine3D_states nstates=$KS local 20 -> 8 A maxits 6, $NP parts x 8 threads"
  simple_exec prg=refine3D_states projfile=p.simple nstates=$KS pose_policy=local lpstart=20 lpstop=8 maxits=6 mskdiam=$MD pgrp=c1 nparts=$NP nthr=8 filt_mode=none > refine.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-160)"; exit 1; }
  RD=$(ls -td $R/disc/*_refine3D_states | head -1); simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $R/refined_fields.txt
  awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' $R/refined_fields.txt > $R/leaves.txt; fi
RD=$(ls -td $R/disc/*_refine3D_states | head -1); log "discovery leaves: $(awk '$1 > 0' $R/leaves.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"
# per-half pose tables of the refined project (half rows -> master rows through meta_rows)
H=$T/simrnasep692
for h in A B; do [ -s $R/poses_$h.txt ] || awk 'NR == FNR {m[FNR] = $1 + 0; n = FNR; next} {l[FNR] = $0} END {for (i = 1; i <= n; i++) print l[m[i]]}' ${H}_$h/meta_rows.txt $R/refined_fields.txt > $R/poses_$h.txt; done
log "per-half poses: $(grep -c . $R/poses_A.txt) / $(grep -c . $R/poses_B.txt) rows"
# ---- (2) leaf tree, (3) merge-up on the real halves
[ -s $R/tree/merge_tree.txt ] || $PY $T/state_merge.py $R/tree --rdir $RD --labels $R/leaves.txt --unfil --link complete --ncut 16,8,4,2 2>&1 | grep -v Warning | grep "^== \|leaf pair" | cut -c1-200
log "merge-up along the leaf tree (real half-stacks)"; [ -s $R/bu/final_units.txt ] || TAG=$RT LEAVES=$R/leaves.txt TREE=$R/tree/merge_tree.txt OUT=$R/bu POSES_A=$R/poses_A.txt POSES_B=$R/poses_B.txt PAR=${PAR:-4} bash $T/bu_merge_real.sbatch 2>&1 | grep -a "units ->\|level\|DONE\|FAILED" | tail -40
U=$R/bu/final_units.txt; M=$(awk '$1 > 0' $U | sort -n | uniq | wc -l); log "$M final states: $(awk '$1 > 0' $U | sort -n | uniq -c | awk '{printf "%d ", $1}')"
# ---- (4) final states -> maps, tree, one cross-half scoring pass, multi-axis sampling
O=$R/coordtree_famall; mkdir -p $O/rec $O/score; cp -f $U $O/units.txt; cd $O/rec
if [ -z "$(ls -d $O/rec/*_reconstruct3D 2>/dev/null)" ]; then cp -f $RD/p.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$U mkdir=no > sel.log 2>&1
  export SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8; simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=$M nparts=$NP nthr=8 > rec.log 2>&1 || { log "final reconstruct3D FAILED"; exit 1; }; fi
RR=$(ls -td $O/rec/*_reconstruct3D | head -1); CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
[ -s $O/final_arm.txt ] || $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/final_arm.txt $RR 2>&1 | grep -v Warning | tail -1
[ -s $O/tree/merge_tree.txt ] || $PY $T/state_merge.py $O/tree --rdir $RR --labels $U --link complete --ncut 2 2>&1 | grep -v Warning | grep "states,\|leaf pair" | cut -c1-160
cd $O/score; if [ -z "$(ls $O/score/*_refine3D/dist_part_neigh_*.dat 2>/dev/null)" ]; then cp -f $O/rec/r.simple s.simple
  XR=$O/score/xref; mkdir -p $XR; VOLS=""
  for s in $(seq 1 $M); do b=recvol_state$(printf %02d $s); cp -f $RR/$b.mrc $XR/$b.mrc; cp -f $RR/${b}_odd.mrc $XR/${b}_even.mrc; cp -f $RR/${b}_even.mrc $XR/${b}_odd.mrc; VOLS="$VOLS vol$s=$XR/$b.mrc"; done
  export SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8; log "scoring all particles against $M maps (cross-half, lp 8, nspace 5000, box 88, $NP parts)"
  simple_exec prg=refine3D projfile=s.simple $VOLS nstates=$M refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=8 maxits=1 center=no nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=$NP nthr=8 > score.log 2>&1 || { log "scoring FAILED: $(grep -a -m2 ERROR score.log | cut -c1-150)"; exit 1; }
  SD=$(ls -td $O/score/*_refine3D | head -1); rm -f $SD/*.mrc; rm -rf $XR; fi
SD=$(ls -td $O/score/*_refine3D | head -1)
# the final states' closed-arm readout keyed by state number (arm_readout lists recvol_stateNN maps)
awk '!/^#/ {n = split($1, p, "state"); s = p[n]; sub(/\.mrc/, "", s); printf "state_%02d.mrc %s %s %s\n", s + 0, $2, $3, $4}' $O/final_arm.txt > $O/final_arm_named.txt
$PY $T/coord_tree2.py $O/nodes2 --tabs "$SD/dist_part_neigh_*.dat" --nstates $M --units $U --order $(seq -s, 1 $M) --tree $O/tree/merge_tree.txt --minp 10000 --arm $O/final_arm_named.txt 2>&1 | grep -v Warning
# ---- (5) the sampled branches' cells: maps + closed-arm readout
RTAG=$RT FAM=all NODESDIR=nodes2 MSKDIAM=$MD bash $T/recon_nodes.sbatch 2>&1 | grep -a "cells\|^   [0-9]\|FAILED"
log DONE

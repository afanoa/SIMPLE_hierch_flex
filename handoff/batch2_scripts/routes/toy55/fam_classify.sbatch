#!/bin/bash
#SBATCH -J klem_famcls
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# within-family K-state classification of a fast run's family: does a multi-state refinement resolve the family's motion, and how far do
# its end states reach? seeds = the bins of one factor (SEED=bins, FAC) or random labels over the family (SEED=random); MODE=states:
# refine3D_states pose_policy=local (the composition stage's recipe); MODE=r3d: refine3D refine=prob_neigh prob_neigh_mode=geom from the
# seeds' reconstructions, optional MASK = pcg_mskfile (hard PCG support: out-of-support density removed). Readout: per state the mean +- sd of
# the true swing (conformers 200-230). env RT, FAM, TAG, SEED (bins), FAC (2), K (8), MODE (states), MASK, LPSTART (15), LPSTOP (8), LP (10), MAXITS (8), LABEL;
# FAMFROM (per-row labels, e.g. comp_k8/refined_states.txt) + FAMLAB ('4+6': the family = rows with these labels; random seeds) instead of a factor
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=$K_/${RT:?}; TAG=${TAG:?}; GT=$T/sim$TAG/gt_particles.txt
SEED=${SEED:-bins}; MODE=${MODE:-states}; K=${K:-8}; FA=$R/factors/fam${FAM:?}_factor${FAC:-2}
O=$R/fcls_f${FAM}_${LABEL:-${MODE}_${SEED}}; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') famcls $RT f$FAM ${LABEL:-$MODE $SEED}] $*"; }
if [ -n "${FAMFROM:-}" ]; then awk -v k=$K -v fl="${FAMLAB:?}" 'BEGIN {srand(7); n = split(fl, a, "+"); for (i = 1; i <= n; i++) in_[a[i]] = 1} {print ($1 in in_) ? 1 + int(k * rand()) : 0}' $FAMFROM > seed_states.txt
elif [ "$SEED" = bins ]; then cp -f $FA/states.txt seed_states.txt
else awk -v k=$K 'BEGIN {srand(7)} {print ($1 > 0) ? 1 + int(k * rand()) : 0}' $FA/states.txt > seed_states.txt; fi
K=$(awk '$1 > 0' seed_states.txt | sort -n | uniq | wc -l)
score(){ awk -v tag="$1" 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; n[s]++; if (k >= 200 && k <= 230) {x = k - 200; sw[s] += x; sq[s] += x * x; ns[s]++}}
  END {m = 0; for (s in n) if (ns[s] > 0) {mu[s] = sw[s] / ns[s]; ord[++m] = s}
    for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (mu[ord[j]] < mu[ord[i]]) {t = ord[i]; ord[i] = ord[j]; ord[j] = t}
    printf "%s (state: mean +- sd swing in A, particles, swing share), by mean:", tag; lo = 99; hi = -1
    for (i = 1; i <= m; i++) {s = ord[i]; sd = sqrt(sq[s] / ns[s] - mu[s] ^ 2); printf " %d:%.1f+-%.1f(%d,%.0f%%)", s, mu[s], sd, n[s], 100 * ns[s] / n[s]
      if (n[s] >= 200 && ns[s] / n[s] >= 0.5) {if (mu[s] < lo) lo = mu[s]; if (mu[s] > hi) hi = mu[s]}}
    printf "\n   %s extent (states >= 200 particles, >= 50%% swing): %.1f .. %.1f A of 0 .. 30\n", tag, lo, hi}' $GT $2; }
score seeds seed_states.txt
cp -f $T/proj$TAG/proj$TAG.simple p.simple
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
if [ "$MODE" = states ]; then
  log "refine3D_states: $K states from the $SEED seeds, local poses, ${LPSTART:-15} -> ${LPSTOP:-8} A, maxits ${MAXITS:-8}"
  simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local pgrp=c1 mskdiam=220 lpstart=${LPSTART:-15} lpstop=${LPSTOP:-8} maxits=${MAXITS:-8} ml_reg=no nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR\|error' ref.log | cut -c1-150)"; exit 1; }
  RD=$(ls -td $O/*_refine3D_states | head -1)
else
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=220 nstates=$K nparts=4 nthr=16 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 'ERROR' rec.log | cut -c1-150)"; exit 1; }
  RR=$(ls -td $O/*_reconstruct3D | head -1); VOLS=""; for s in $(seq 1 $K); do v=$(ls $RR/recvol_state$(printf %02d $s).mrc 2>/dev/null); [ -n "$v" ] || { log "no seed map for state $s in $RR"; exit 1; }; VOLS="$VOLS vol$s=$v"; done
  MK=""; [ -n "${MASK:-}" ] && MK="rec_backend=pcg pcg_mskfile=$MASK"
  log "refine3D prob_neigh geom: $K states from the $SEED seeds, lp ${LP:-10}, maxits ${MAXITS:-8} $MK"
  simple_exec prg=refine3D projfile=p.simple $VOLS nstates=$K refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=220 lp=${LP:-10} maxits=${MAXITS:-8} center=no $MK nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 'ERROR' ref.log | cut -c1-150)"; exit 1; }
  RD=$(ls -td $O/*_refine3D | head -1)
fi
PF=$RD/p.simple; [ -s $PF ] || PF=$O/p.simple
simple_exec prg=print_project_field projfile=$PF oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/final_states.txt
score final final_states.txt
log DONE

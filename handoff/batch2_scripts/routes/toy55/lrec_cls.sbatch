#!/bin/bash
#SBATCH -J lrec_cls
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=120G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/logs/%x.%j.out
# leaf-discovery recipe tests on an oracle family of the simswing30 benchmark: K-state refine3D_states on the swing family (or a
# random subset of it), seeded at random, by true-swing bins (ceiling) or from a per-row label file; readout per state mean +- sd of
# the true swing (conformers 200-230) and eta2 (share of the family's swing variance explained by the states).
# env LABEL (run name), K (8), SEED (random|truth|file), SEEDFILE, NSUB (random subset of the family; all), TAG (swing30),
# FAMFROM (per-row family labels) FAMLAB (1), LPSTART (15), LPSTOP (8), MAXITS (8), EXTRA (more refine3D_states keys), SUBSEED (7)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; KD=$A/rnasep_klem; TAG=${TAG:-swing30}; GT=$T/sim$TAG/gt_particles.txt
FAMFROM=${FAMFROM:-$KD/ff2_swing30/oracle_swing_rows.txt}; FAMLAB=${FAMLAB:-1}; K=${K:-8}; SEED=${SEED:-random}
O=$KD/lrec30/${LABEL:?}; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') lrec $LABEL] $*"; }
# family rows (optionally a random subset of NSUB), then the seed labels over them
awk -v fl="$FAMLAB" -v ns="${NSUB:-0}" -v sd="${SUBSEED:-7}" 'BEGIN {n = split(fl, a, "+"); for (i = 1; i <= n; i++) in_[a[i]] = 1}
  NR == FNR {if ($1 in in_) nf++; next} FNR == 1 {srand(sd); p = (ns > 0) ? ns / nf : 1}
  {print (($1 in in_) && rand() < p) ? 1 : 0}' $FAMFROM $FAMFROM > fam_rows.txt
case $SEED in
  random) awk -v k=$K 'BEGIN {srand(11)} {print ($1 > 0) ? 1 + int(k * rand()) : 0}' fam_rows.txt > seed_states.txt;;
  truth)  awk -v k=$K 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; x = c[r] - 200; print ($1 > 0 && x >= 0 && x <= 30) ? 1 + int(x * k / 31) : 0}' $GT fam_rows.txt > seed_states.txt;;
  file)   paste fam_rows.txt ${SEEDFILE:?} | awk '{print ($1 > 0) ? $2 : 0}' > seed_states.txt;;
  *) log "unknown SEED $SEED"; exit 1;;
esac
K=$(awk '$1 > 0' seed_states.txt | sort -n | uniq | wc -l)
score(){ awk -v tag="$1" 'NR == FNR {if ($0 !~ /^#/) c[$1] = $2; next} {r++; s = $1; if (s <= 0) next; k = c[r]; n[s]++; if (k >= 200 && k <= 230) {x = k - 200; sw[s] += x; sq[s] += x * x; ns[s]++; T1 += x; T2 += x * x; NT++}}
  END {m = 0; for (s in n) if (ns[s] > 0) {mu[s] = sw[s] / ns[s]; ord[++m] = s; W += sq[s] - sw[s] ^ 2 / ns[s]}
    for (i = 1; i <= m; i++) for (j = i + 1; j <= m; j++) if (mu[ord[j]] < mu[ord[i]]) {t = ord[i]; ord[i] = ord[j]; ord[j] = t}
    printf "%s: %d swing particles, eta2 %.3f; per state mean +- sd swing (particles), by mean:", tag, NT, 1 - W / (T2 - T1 ^ 2 / NT); lo = 99; hi = -1
    for (i = 1; i <= m; i++) {s = ord[i]; sd = sqrt(sq[s] / ns[s] - mu[s] ^ 2); printf " %d:%.1f+-%.1f(%d)", s, mu[s], sd, n[s]
      if (n[s] >= 200 && ns[s] / n[s] >= 0.5) {if (mu[s] < lo) lo = mu[s]; if (mu[s] > hi) hi = mu[s]}}
    printf "\n   %s extent: %.1f .. %.1f A of 0 .. 30\n", tag, lo, hi}' $GT $2; }
score seeds seed_states.txt
cp -f $T/proj$TAG/proj$TAG.simple p.simple
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/seed_states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1
log "refine3D_states: $K states, $SEED seeds, $(awk '$1 > 0' seed_states.txt | wc -l) particles, local poses, ${LPSTART:-15} -> ${LPSTOP:-8} A, maxits ${MAXITS:-8} ${EXTRA:-}"
simple_exec prg=refine3D_states projfile=p.simple nstates=$K pose_policy=local pgrp=c1 mskdiam=220 lpstart=${LPSTART:-15} lpstop=${LPSTOP:-8} maxits=${MAXITS:-8} ml_reg=no nparts=4 nthr=16 ${EXTRA:-} > ref.log 2>&1 || { log "refine3D_states FAILED: $(grep -a -m2 'ERROR\|error' ref.log | cut -c1-150)"; exit 1; }
RD=$(ls -td $O/*_refine3D_states | head -1)
simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/final_states.txt
grep -a "STATE OVERLAP" ref.log | awk '{printf " %s", $NF} END {print ""}' | sed 's/^/state overlap per iteration:/'
score final final_states.txt
log DONE

#!/bin/bash
#SBATCH -J klem_smdiag10
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=160G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# split-test diagnostic, round 10 = splitbench16 (pairs tested at the equal N of their smaller member); round 8 = splitbench15 (docking to the other half's consensus, merged groups re-docked, SSNR junk); round 7 = splitbench14 (merged groups rebuilt at N_eq) --dock 2; round 6 = splitbench13 --dock 2 (test maps docked on the tested band; EM references on the EM band); m4 = re-test
# band 4 shells beyond the EM band. Round 4 (splitbench12 --dock: group maps translated in 3-D to their half consensus in the EM references and the test):
# every variant = cross-half EM (xem 3), re-test beyond the EM band (hold-out), stratified null, per-shell log, junk screen; consensus poses
# unless 'comb' (each leaf's state-refined poses); fulln (groups at full size) or fix (equal N).
# env CASES, VARS, PAR (3), JP (0.75), JS (0.35); out rnasep_klem/smdiag/<case>/<variant>/{bench.log,summary.txt}
set -uo pipefail
export OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
O=$K_/smdiag; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') smdiag10] $*"; }
[ -s $O/cases.txt ] || { log "no $O/cases.txt (written by smdiag.sbatch)"; exit 1; }
CASES=${CASES:-null50 two50 r11 r12 r21 r22}; VARS=${VARS:-g16_fix_j}
run(){ cs=$1; v=$2; set -- $(grep "^$cs " $O/cases.txt); S=$2; PC=$3; PB=$4; U=$5; E=$6; CP=$7; GT=$8; d=$O/$cs/$v; mkdir -p $d
  [ -s $d/summary.txt ] && return 0
  case $v in *comb*) [ "$PB" = - ] && return 0; PC=$PB;; esac
  X="--xem 3 --retest_beyond_em 1 --em_min 1000 --shells 1 --dock 2"; case $v in *_m4*) X="$X --rb_margin 4";; esac; case $v in *fulln*) X="$X --fulln 1";; esac; case $v in *sm3*) X="$X --specmatch 3";; esac
  case $v in *_j) X="$X --junk_proj ${JP:-0.75} --junk_ssnr ${JS:-0.35}";; esac
  [ -s $O/$cs/cand_B.txt ] || $PY -c "
import numpy as np; eo = np.loadtxt('$E', dtype=int); u = np.loadtxt('$U', dtype=int)
np.savetxt('$O/$cs/cand_A.txt', np.where(eo == 0, u, 0), fmt='%d'); np.savetxt('$O/$cs/cand_B.txt', np.where(eo == 1, u, 0), fmt='%d')"
  t0=$(date +%s); $PY -u $T/splitbench16.py --sim $S,$S --poses $PC,$PC --cand $O/$cs/cand_A.txt,$O/$cs/cand_B.txt --out $d --cache $CP --ctfabs 0 --pair identity --softem 0 --arch 0 --strata 1 \
    --gt $GT,$GT --summary $d/summary.txt $X > $d/bench.log 2>&1 && echo "   $cs $v ($(( $(date +%s) - t0 )) s): $(head -1 $d/summary.txt)" || echo "   $cs $v FAILED: $(grep -a -m1 -A3 'Error\|Traceback' $d/bench.log | tail -2)"; }
log "cases $CASES; variants $VARS"
for cs in $CASES; do mkdir -p $O/$cs; for v in $VARS; do
  while [ $(jobs -rp | wc -l) -ge ${PAR:-4} ]; do sleep 5; done; run $cs $v &
  [ ! -s $O/$cs/cand_B.txt ] && wait
done; done; wait
log "SUMMARY (K* sizes NMI)"; for cs in $CASES; do for v in $VARS; do [ -s $O/$cs/$v/summary.txt ] && echo "$cs $v: $(head -1 $O/$cs/$v/summary.txt)"; done; done
log DONE

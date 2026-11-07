#!/bin/bash
#SBATCH -J klem_smdiag12
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=160G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# split-test diagnostic, round 12: power for low-resolution / localized differences (the real-data closed-arm state 6 of rflex_f692u7 node_r
# merged by the 16.7-10.2 A re-test). splitbench16 variants: base (EM <= 20 A, re-test from 16.7 A), em30 (EM <= 30 A, re-test from 23 A),
# foc (comparison inside the node's focus mask, re-test 4 shells beyond the EM band), em30foc (both). Real f692u7 (read only: half-stacks,
# node poses, candidate files, rsmask) and controls rflex_s30a r11 (must merge) / r22 (must keep) at their own refined poses and masks.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K_/smdiag; log(){ echo "[$(date '+%F %T') smdiag12] $*"; }
U7=$K_/rflex_f692u7/node_r; H=$T/simrnasep692; R=$K_/rflex_s30a
var(){ case $1 in base) echo "";; em30) echo "--em_lim 30";; foc) echo "--focus $2 --rb_margin 4";; em30foc) echo "--em_lim 30 --focus $2 --rb_margin 4";; esac; }
run(){ cs=$1; v=$2; d=$O/$cs/$v; mkdir -p $d; [ -s $d/summary.txt ] && return 0; t0=$(date +%s)
  case $cs in
    f692u7) $PY -u $T/splitbench16.py --sim ${H}_A,${H}_B --poses $U7/poses_A.txt,$U7/poses_B.txt --cand $U7/ftest/cand_A.txt,$U7/ftest/cand_B.txt --out $d --cache $T/benchcache_rnasep692 --ctfabs 1 \
              --pair identity --softem 0 --arch 0 --strata 1 --xem 3 --retest_beyond_em 1 --em_min 1000 --dock 2 --junk_proj 0.75 --junk_ssnr 0.35 --shells 1 $(var $v $U7/rsmask.mrc) --summary $d/summary.txt > $d/bench.log 2>&1;;
    r11|r22) S=$T/simswing30; $PY -u $T/splitbench16.py --sim $S,$S --poses $R/post/poses_combined.txt,$R/post/poses_combined.txt --cand $O/$cs/cand_A.txt,$O/$cs/cand_B.txt --out $d --cache $T/benchcache_full_swing30 --ctfabs 0 \
              --pair identity --softem 0 --arch 0 --strata 1 --xem 3 --retest_beyond_em 1 --em_min 1000 --dock 2 --junk_proj 0.75 --junk_ssnr 0.35 --shells 1 --gt $S/gt_particles.txt,$S/gt_particles.txt $(var $v $R/node_$cs/mask.mrc) --summary $d/summary.txt > $d/bench.log 2>&1;;
  esac && echo "   $cs $v ($(( $(date +%s) - t0 )) s): $(head -1 $d/summary.txt | cut -c1-200)" || echo "   $cs $v FAILED: $(grep -a -m1 -A3 'Error\|Traceback' $d/bench.log | tail -2)"
  if [ $cs = f692u7 ] && [ -s $d/summary.txt ]; then $PY - <<PYEOF
import numpy as np
ua, ub = np.loadtxt('$U7/ftest/cand_A.txt', dtype=int), np.loadtxt('$U7/ftest/cand_B.txt', dtype=int); ca, cb = np.loadtxt('$d/children_A.txt', dtype=int), np.loadtxt('$d/children_B.txt', dtype=int)
u, c = np.r_[ua, ub], np.r_[ca, cb]; to = {}
for s in np.unique(u[u > 0]): x = c[u == s]; x = x[x > 0]; to[int(s)] = int(np.bincount(x).argmax()) if len(x) else 0
grp = {}
for s, g in to.items(): grp.setdefault(g, []).append(s)
print('      $cs $v: units -> supported groups (plurality): %s; closed unit 6 %s' % ('  '.join('%d:{%s}' % (g, ','.join(map(str, m))) for g, m in sorted(grp.items())), 'ALONE' if len(grp.get(to.get(6, -1), [])) == 1 else 'merged with %s' % grp.get(to.get(6, -1))))
PYEOF
  fi; }
log "start"
for cs in f692u7 r11 r22; do for v in base em30 foc em30foc; do
  [ $cs != f692u7 ] && [ $v = base ] && continue
  while [ $(jobs -rp | wc -l) -ge 4 ]; do sleep 5; done; run $cs $v &
done; done; wait
log DONE

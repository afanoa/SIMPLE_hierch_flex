#!/bin/bash
#SBATCH -J klem_smdiag
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=160G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# split-test diagnostic (splitbench10): per node, its candidate states tested on the even/odd halves of the one stack with the cross-half EM
# (xem 3, stratified null); variants: poses (cons = the input consensus poses, comb = each leaf's own state-refined poses), rb = re-test only
# beyond the EM band (frequency hold-out), sm1 = per-shell spectrum matching. Cases: rflex_s30a nodes (truth simswing30) and the null panels.
# env CASES (default all), VARS (default all), PAR (3); out rnasep_klem/smdiag/<case>/<variant>/{bench.log,summary.txt}
set -uo pipefail
export OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
O=$K_/smdiag; R=$K_/rflex_s30a; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') smdiag] $*"; }
# case: name sim poses_cons poses_comb units eo cache gt
cat > $O/cases.txt <<EOF
null50 $T/simnull_50k $T/lat_null_50k/ptcl3D_r0.txt - $K_/panel_null_50k/refined_states.txt $K_/panel_null_50k/bu/eo.txt $T/benchcache_full_null_50k $T/simnull_50k/gt_particles.txt
two50 $T/simtwo_50k $T/lat_two_50k/ptcl3D_r0.txt - $K_/panel_two_50k/refined_states.txt $K_/panel_two_50k/bu/eo.txt $T/benchcache_full_two_50k $T/simtwo_50k/gt_particles.txt
r11 $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt $R/post/poses_combined.txt $R/node_r11/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
r12 $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt $R/post/poses_combined.txt $R/node_r12/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
r21 $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt $R/post/poses_combined.txt $R/node_r21/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
r22 $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt $R/post/poses_combined.txt $R/node_r22/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
r1 $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt - $R/node_r1/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
r2 $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt - $R/node_r2/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
r $T/simswing30 $T/lat_swing30/ptcl3D_r0.txt - $R/node_r/refined_states.txt $R/post/bu/eo.txt $T/benchcache_full_swing30 $T/simswing30/gt_particles.txt
EOF
CASES=${CASES:-null50 two50 r11 r12 r21 r22 r1 r2 r}; VARS=${VARS:-cons cons_rb cons_sm1 cons_rb_sm1 comb_rb_sm1}
run(){ cs=$1; v=$2; set -- $(grep "^$cs " $O/cases.txt); S=$2; PC=$3; PB=$4; U=$5; E=$6; CP=$7; GT=$8; d=$O/$cs/$v; mkdir -p $d
  [ -s $d/summary.txt ] && return 0
  P=$PC; case $v in comb*) [ "$PB" = - ] && return 0; P=$PB;; esac
  X="--xem 3"; case $v in *_rb*) X="$X --retest_beyond_em 1";; esac; case $v in *sm1*) X="$X --specmatch 1";; *sm2*) X="$X --specmatch 2";; esac
  [ -s $O/$cs/cand_B.txt ] || $PY -c "
import numpy as np; eo = np.loadtxt('$E', dtype=int); u = np.loadtxt('$U', dtype=int)
np.savetxt('$O/$cs/cand_A.txt', np.where(eo == 0, u, 0), fmt='%d'); np.savetxt('$O/$cs/cand_B.txt', np.where(eo == 1, u, 0), fmt='%d')"
  t0=$(date +%s); $PY -u $T/splitbench10.py --sim $S,$S --poses $P,$P --cand $O/$cs/cand_A.txt,$O/$cs/cand_B.txt --out $d --cache $CP --ctfabs 0 --pair identity --softem 0 --arch 0 --strata 1 \
    --gt $GT,$GT --summary $d/summary.txt $X > $d/bench.log 2>&1 && echo "   $cs $v ($(( $(date +%s) - t0 )) s): $(head -1 $d/summary.txt)" || echo "   $cs $v FAILED: $(grep -a -m1 -A3 'Error\|Traceback' $d/bench.log | tail -2)"; }
log "cases $CASES; variants $VARS"
for cs in $CASES; do mkdir -p $O/$cs; for v in $VARS; do
  case $cs in r|r1|r2) case $v in cons_rb|cons_rb_sm1) ;; *) continue;; esac;; esac   # large nodes: the hold-out variants only
  while [ $(jobs -rp | wc -l) -ge ${PAR:-3} ]; do sleep 5; done; run $cs $v &
  [ ! -s $O/$cs/cand_B.txt ] && wait   # first variant of a case writes the candidate files
done; done; wait
log "SUMMARY (K* sizes NMI)"; for cs in $CASES; do for v in $VARS; do [ -s $O/$cs/$v/summary.txt ] && echo "$cs $v: $(head -1 $O/$cs/$v/summary.txt)"; done; done
log DONE

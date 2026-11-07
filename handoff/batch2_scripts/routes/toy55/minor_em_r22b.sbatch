#!/bin/bash
#SBATCH -J klem_minor_em_b
#SBATCH -p norm
#SBATCH -c 16
#SBATCH --mem=72G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# minority states on real r22, the general form: the node's refined states (12, populations evened by hard refinement) PLUS the two end bins
# of the equal-width arm-axis seed as extra components 13 / 14 (their particles start there), soft vs hard assignment at fixed poses
# (minor_em.py, arm-mask deviations); closed-arm readout of every state map at the start and after the EM. Inputs from minor_em_r22.sbatch.
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
FA=$K_/focused_arm_f692rf_r22; O=$K_/minor_em_r22; cd $O; log(){ echo "[$(date '+%F %T') minor-em r22b] $*"; }
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc
for h in A B; do [ -s $O/plusends_$h.txt ] || paste $O/refined_$h.txt $O/seed_$h.txt | awk '{print ($1 == 0) ? 0 : ($2 == 1 ? 13 : ($2 == 12 ? 14 : $1))}' > $O/plusends_$h.txt; done
log "start labels $(cat $O/plusends_A.txt $O/plusends_B.txt | awk '$1 > 0' | sort -n | uniq -c | awk '{printf "%d ", $1}')"
run(){ tag=$1; soft=$2; d=$O/$tag; mkdir -p $d; t0=$(date +%s)
  $PY -u $T/minor_em.py --sim ${H}_A,${H}_B --poses $O/poses_A.txt,$O/poses_B.txt --labels $O/plusends_A.txt,$O/plusends_B.txt --cache $T/benchcache_rnasep692 \
      --ctfabs 1 --cons_ref $CONS --mask $FA/arm_mask.mrc --soft $soft --iters ${ITERS:-8} --em_lim ${EMLIM:-12} --out $d > $d/minor_em.log 2>&1
  log "$tag: exit $? in $(( $(date +%s) - t0 )) s"
  for k in start soft hard; do [ -d $d/maps_$k ] && $PY $T/arm_readout.py $CONS $CL $d/arm_$k.txt $d/maps_$k 2>&1 | grep -v Warning | tail -1 | sed "s/^/   $tag $k: /"; done; }
run plusends_soft 1 & sleep 60; run plusends_hard 0 & wait
log "SUMMARY (closed-arm/open-arm per state map, by closed-arm; occupancies from the last iteration line)"
for r in plusends_soft plusends_hard; do d=$O/$r; echo "== $r"; grep -a '^start\|^iteration 8' $d/minor_em.log | cut -c1-400
  for k in start soft hard; do [ -s $d/arm_$k.txt ] && echo "   $k: $(grep -v '^#' $d/arm_$k.txt | sort -k2,2gr | awk '{n = split($1, p, "/"); printf "%s %.2f/%.2f  ", p[n], $2, $3}')"; done; done
log DONE

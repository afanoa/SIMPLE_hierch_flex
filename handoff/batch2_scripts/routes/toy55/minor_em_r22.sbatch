#!/bin/bash
#SBATCH -J klem_minor_em
#SBATCH -p norm
#SBATCH -c 24
#SBATCH --mem=120G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# minority states on real r22 (f692rf): soft vs hard assignment at fixed (parent) poses, cross-half references, arm-mask deviations, per-shell
# reference scale (minor_em.py), started from (a) the equal-width arm-axis bins (children_w1: small end bins, the 12-bin end read 0.81) and
# (b) the arm-focused refine3D states seeded by them (populations evened to 27-35k, closed-arm max 0.63); closed-arm readout of every state map
# at the start and after the EM. out rnasep_klem/minor_em_r22/<run>/{minor_em.log,maps_*,arm_*.txt}
set -uo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; H=$T/simrnasep692
CK=$K_/cheapkids_f692rf_r22; FA=$K_/focused_arm_f692rf_r22; O=$K_/minor_em_r22; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') minor-em r22] $*"; }
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc
for h in A B; do mr=${H}_$h/meta_rows.txt
  [ -s $O/poses_$h.txt ] || awk 'NR == FNR {l[FNR] = $0; next} {print l[$1 + 0]}' $CK/ptcl3D_parent.txt $mr > $O/poses_$h.txt
  [ -s $O/seed_$h.txt ] || awk 'NR == FNR {l[FNR] = $1; next} {print l[$1 + 0] + 0}' $CK/children_w1.txt $mr > $O/seed_$h.txt
  [ -s $O/refined_$h.txt ] || awk 'NR == FNR {l[FNR] = $1; next} {print l[$1 + 0] + 0}' $FA/states.txt $mr > $O/refined_$h.txt
done
log "half rows: A $(wc -l < $O/poses_A.txt) B $(wc -l < $O/poses_B.txt); seed labels $(cat $O/seed_A.txt $O/seed_B.txt | awk '$1 > 0' | sort -n | uniq -c | awk '{printf "%d ", $1}')"
run(){ tag=$1; lab=$2; soft=$3; d=$O/$tag; mkdir -p $d; t0=$(date +%s)
  $PY -u $T/minor_em.py --sim ${H}_A,${H}_B --poses $O/poses_A.txt,$O/poses_B.txt --labels $O/${lab}_A.txt,$O/${lab}_B.txt --cache $T/benchcache_rnasep692 \
      --ctfabs 1 --cons_ref $CONS --mask $FA/arm_mask.mrc --soft $soft --iters ${ITERS:-8} --em_lim ${EMLIM:-12} --out $d > $d/minor_em.log 2>&1
  log "$tag: exit $? in $(( $(date +%s) - t0 )) s"
  for k in start soft hard; do [ -d $d/maps_$k ] && $PY $T/arm_readout.py $CONS $CL $d/arm_$k.txt $d/maps_$k 2>&1 | grep -v Warning | tail -1 | sed "s/^/   $tag $k: /"; done; }
run seed_soft seed 1 & sleep 60; run seed_hard seed 0 & sleep 60; run refined_soft refined 1 & wait
log "SUMMARY (closed-arm, open-arm per state map, by closed-arm; occupancies from the last iteration line)"
for r in seed_soft seed_hard refined_soft; do d=$O/$r; echo "== $r: $(grep -a 'frame:' $d/minor_em.log | cut -c1-120)"; grep -a '^start\|^iteration 8' $d/minor_em.log | cut -c1-400
  for k in start soft hard; do [ -s $d/arm_$k.txt ] && echo "   $k: $(grep -v '^#' $d/arm_$k.txt | sort -k2,2gr | awk '{n = split($1, p, "/"); printf "%s %.2f/%.2f  ", p[n], $2, $3}')"; done; done
log DONE

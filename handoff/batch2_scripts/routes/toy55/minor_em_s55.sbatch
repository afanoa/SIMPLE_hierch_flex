#!/bin/bash
#SBATCH -J klem_minor_s55
#SBATCH -p norm
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# minority states with truth (sim55_200k: 9.3 % closed over 7 conformers, 20 % junk): 12 equal-width bins of a noisy closed score (d' 1.47 =
# AUC 0.85, the real per-particle separability) as the seed, then soft vs hard assignment at fixed poses (minor_em.py), with and without the
# r22 arm mask on the deviations; readout = closed share per state, purity / recall of the most closed state, per iteration.
set -uo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; S=$T/sim55_200k
O=$K_/minor_em_s55; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') minor-em s55] $*"; }; DP=${DPRIME:-1.47}
[ -s $O/seed_B.txt ] || $PY $T/noisy_seeds.py ${S}_A/gt_particles.txt ${S}_B/gt_particles.txt $O/seed_A.txt $O/seed_B.txt --dprime $DP --nbins 12
run(){ tag=$1; soft=$2; mk=$3; d=$O/$tag; mkdir -p $d; t0=$(date +%s); MX=""; [ $mk = 1 ] && MX="--mask $K_/focused_arm_f692rf_r22/arm_mask.mrc"
  $PY -u $T/minor_em.py --sim ${S}_A,${S}_B --poses $T/lat_55_200k_A/ptcl3D_r0.txt,$T/lat_55_200k_B/ptcl3D_r0.txt --labels $O/seed_A.txt,$O/seed_B.txt \
      --cache $T/benchcache_55_200k --ctfabs 0 --cons_ref $T/maps160/toy_consensus.mrc $MX --soft $soft --iters ${ITERS:-8} --em_lim ${EMLIM:-12} \
      --gt ${S}_A/gt_particles.txt,${S}_B/gt_particles.txt --out $d > $d/minor_em.log 2>&1
  log "$tag: exit $? in $(( $(date +%s) - t0 )) s; $(grep -a '^iteration 8' $d/minor_em.log | grep -o 'most closed state.*')"; }
run soft_mask 1 1 & run hard_mask 0 1 & run soft_nomask 1 0 & run hard_nomask 0 0 & wait
log "SUMMARY"; for r in soft_mask hard_mask soft_nomask hard_nomask; do echo "== $r: $(grep -a 'frame:' $O/$r/minor_em.log | cut -c1-110)"; grep -a '^start\|^iteration' $O/$r/minor_em.log | grep -o '^[a-z]* *[0-9]*\|closed share.*' | paste - - | cut -c1-300; done
log DONE

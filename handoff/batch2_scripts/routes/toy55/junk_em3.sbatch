#!/bin/bash
#SBATCH -J klem_junk_em3
#SBATCH -p norm
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# junk component prototype v3 (junk_em3.py: klem start maps, free maps reconstructed like the states, per-shell reference scale) on the synthetic junk sets (null_50k: one conformer + 10k junk; two_50k: two conformers + 10k junk),
# fixed refine3D poses and states. Configurations: noise component only / one free component only / both; hide 1 = the truth junk starts
# inside random real states, hide 0 = the refined states as they are. env SETS (default null_50k two_50k), CFGS, PAR (6)
set -uo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem
O=$K_/junk_em3; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') junk-em] $*"; }
SETS=${SETS:-null_50k two_50k}; CFGS=${CFGS:-k3_h1 k3_h0 k3nb_h1 low3_h1}
run(){ st=$1; cf=$2; d=$O/${st}_$cf; mkdir -p $d; [ -s $d/calls_B.txt ] && return 0
  case $cf in k3_*) f="--free 2 --free_init klem --beta 1";; k3nb_*) f="--free 2 --free_init klem --beta 0";; low3_*) f="--free 1 --free_init low --beta 1";; esac
  h=${cf##*_h}; t0=$(date +%s)
  $PY -u $T/junk_em3.py --sim $T/sim$st --poses $T/lat_$st/ptcl3D_r0.txt --eo $K_/panel_$st/bu/eo.txt --labels $K_/panel_$st/refined_states.txt \
      --gt $T/sim$st/gt_particles.txt --cache $T/benchcache_full_$st --out $d $f --hide $h --iters 8 > $d/junk_em3.log 2>&1
  log "$st $cf: exit $? in $(( $(date +%s) - t0 )) s; $(grep 'iteration 8' $d/junk_em3.log | cut -c1-700)"; }
log "sets $SETS; configurations $CFGS"
for st in $SETS; do for cf in $CFGS; do
  while [ $(jobs -rp | wc -l) -ge ${PAR:-8} ]; do sleep 5; done; run $st $cf &
done; done; wait
log "SUMMARY"; for st in $SETS; do for cf in $CFGS; do d=$O/${st}_$cf; echo "== $st $cf"; grep 'hidden\|lowest-scoring\|junk types\|start |\|iteration 8\|reference scale\|final components' $d/junk_em3.log | cut -c1-500; done; done
log DONE

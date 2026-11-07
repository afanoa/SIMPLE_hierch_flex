#!/bin/bash
#SBATCH -J klem_fastfac
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# fast factor test on a synthetic set (one node, no worker jobs, no refinement: the simulated poses are near-true): NSUB particles,
# one whole-particle board, regions inside the periphery (derive_regions2 --within), one focused board per region (4 at a time),
# then factor_sample.py (dependence between the region boards, factors, even bins) with the truth readout.
# env TAG (simulation: sim<TAG>, proj<TAG>, lat_<TAG>), RT (ff_<TAG>), NSUB (100000) or SUBSETFILE (rows), KREG (8), CODES (truth codes), NBIN (8)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; W=$A/routes/toy55_f692v2; K=$A/rnasep_klem
TAG=${TAG:?}; R=$K/${RT:-ff_$TAG}; D=$R/fam1; mkdir -p $R/expand $D/refine/1_refine3D_states $D/v2; cd $R; log(){ echo "[$(date '+%F %T') fastfac $TAG] $*"; }
export FB=$K/f692e/klem2_bin CTFABS=0 PSCALE=0 SERVER_LIST= NPARTS=
NR=$(grep -c . $T/lat_$TAG/ptcl3D_r0.txt); [ -n "${SUBSETFILE:-}" ] && [ ! -s $R/expand/family1_subset.txt ] && cp -f $SUBSETFILE $R/expand/family1_subset.txt   # e.g. an oracle composition class
[ -s $R/expand/family1_subset.txt ] || $PY -c "import numpy as np; r = np.sort(np.random.default_rng(1).choice($NR, min($NR, ${NSUB:-100000}), replace=False)) + 1; np.savetxt('$R/expand/family1_subset.txt', r, fmt='%d')"
[ -s $D/meta_fb.npz ] || $PY $T/make_toy_meta.py $T/lat_$TAG/ptcl3D_r0.txt $D/meta_fb.npz > /dev/null
[ -s $D/refine/1_refine3D_states/p.simple ] || cp -f $T/proj$TAG/proj$TAG.simple $D/refine/1_refine3D_states/p.simple
PS=$D/refine/1_refine3D_states/p.simple; SUB=$R/expand/family1_subset.txt; log "$(grep -c . $SUB) of $NR particles"
fb(){ local out=$1 mk=$2 lps=$3 sig=$4 nthr=$5; [ -s $out/klem_latent.txt ] && return 0
  NTHR=$nthr N1=7 N2=7 NFREE=1 FOCUSMASK=$mk META=$D/meta_fb.npz bash $W/fboard.sh $out $PS $SUB "$lps" "$sig" > $out.out 2>&1 || log "board FAILED $out: $(tail -1 $out.out | cut -c1-120)"; }
# ---- whole-particle board, periphery, regions
fb $D/klem_fb "" "14 12 10 8 8" "1.2 1.0 0.8 0.7 0.7" 64; log "whole-particle board done"
[ -s $D/v2/periphery.mrc ] || $PY $T/periphery_mask.py $T/maps160/toy_consensus.mrc $D/v2/periphery.mrc --core ${CORE:-0.40} 2>&1 | grep -v Warning | tail -1
[ -s $D/v2/cand/region_1.mrc ] || $PY $T/derive_regions2.py $D/klem_fb $D/v2/cand --k ${KREG:-8} --maxreg ${KREG:-8} --within $D/v2/periphery.mrc 2>&1 | grep -v Warning | tail -3
# ---- one focused board per region (4 at a time x 16 threads)
for f in $D/v2/cand/region_*.mrc; do r=$(basename $f .mrc); r=${r#region_}; mkdir -p $D/region_c$r
  while [ $(jobs -rp | wc -l) -ge 4 ]; do sleep 5; done; fb $D/region_c$r/board $f "12 10 10 8 8" "1.2 1.0 1.0 0.8 0.8" 16 &
done; wait; log "$(ls -d $D/region_c*/board/klem_latent.txt 2>/dev/null | wc -l) region boards done"
# ---- factors with the truth readout
$PY $T/factor_sample.py $R $R/factors --fams 1 --glob 'region_c*' --nbin ${NBIN:-8} --budget ${BUDGET:-1500} --minend ${MINEND:-500} --gt $T/sim$TAG/gt_particles.txt --codes ${CODES:?} 2>&1 | grep -v Warning
log DONE

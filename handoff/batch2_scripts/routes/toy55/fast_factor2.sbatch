#!/bin/bash
#SBATCH -J klem_fastfac2
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# fast factor test with COMPOSITION FIRST (one node, no worker jobs, near-true simulated poses): NSUB particles -> whole-particle
# board -> families = spectral clusters of its cells (number from the largest eigengap, NFAM to force; every on-grid particle -> its
# cell's family) -> regions inside the periphery from the whole board (spatial, derived once) -> per family one focused board per region ->
# factor_sample.py with the truth readout and org_readout.py (top split vs composition).
# env TAG, RT (ff2_<TAG>), NSUB (100000), NFAM (eigengap), FAMFROM (per-row family labels, e.g. comp_k8/refined_states.txt: skip the split), KREG (6), NFREEB (0: on the synthetic sets the free component swallowed
# 57-64 % of the real particles and 1 % of the junk), CODES (truth codes, passed through the environment), NBIN (8)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; W=$A/routes/toy55_f692v2; K=$A/rnasep_klem
TAG=${TAG:?}; R=$K/${RT:-ff2_$TAG}; mkdir -p $R/expand $R/whole; cd $R; log(){ echo "[$(date '+%F %T') fastfac2 $TAG] $*"; }
export FB=$K/f692e/klem2_bin CTFABS=0 PSCALE=0 SERVER_LIST= NPARTS=
NR=$(grep -c . $T/lat_$TAG/ptcl3D_r0.txt); ALL=$R/expand/all_subset.txt
[ -s $ALL ] || $PY -c "import numpy as np; r = np.sort(np.random.default_rng(1).choice($NR, min($NR, ${NSUB:-100000}), replace=False)) + 1; np.savetxt('$ALL', r, fmt='%d')"
META=$R/meta_toy.npz; [ -s $META ] || $PY $T/make_toy_meta.py $T/lat_$TAG/ptcl3D_r0.txt $META > /dev/null
PS=$R/whole/p.simple; [ -s $PS ] || cp -f $T/proj$TAG/proj$TAG.simple $PS; log "$(grep -c . $ALL) of $NR particles"
fb(){ local out=$1 mk=$2 sub=$3 lps=$4 sig=$5 nthr=$6; [ -s $out/klem_latent.txt ] && return 0
  NTHR=$nthr N1=7 N2=7 NFREE=${NFREEB:-0} FOCUSMASK=$mk META=$META bash $W/fboard.sh $out $PS $sub "$lps" "$sig" > $out.out 2>&1 || log "board FAILED $out: $(tail -1 $out.out | cut -c1-120)"; }
# ---- (1) composition: whole-particle board, families from its cell graph
fb $R/whole/board "" $ALL "20 14 12 10 10" "1.5 1.2 1.0 1.0 1.0" 64; log "whole-particle board done"
if [ -n "${FAMFROM:-}" ]; then $PY -c "import numpy as np; L = np.loadtxt('$FAMFROM', dtype=int); k = int(L.max()); [np.savetxt('$R/expand/family%d_subset.txt' % j, np.where(L == j)[0] + 1, fmt='%d') for j in range(1, k + 1)]; open('$R/expand/nfam.txt', 'w').write('%d\\n' % k); print('families from %s: %d' % ('$FAMFROM', k))"; else
$PY - <<PYEOF 2>&1 | grep -v Warning   # families = spectral clusters of the whole board's cells (number from the largest eigengap unless NFAM); every on-grid particle -> its cell's family
import glob, numpy as np
def rd(p): h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23])).reshape(n, n, n).astype(float)
fs = sorted(glob.glob('$R/whole/board/cells/cell_*.mrc')); V = np.array([rd(f) for f in fs]); D = V - V.mean(0); X = np.array([np.fft.fftn(d).ravel() for d in D]); X /= np.linalg.norm(X, axis=1)[:, None]
S = np.maximum(np.real(X @ X.conj().T), 0); np.fill_diagonal(S, 0); dg = S.sum(1) + 1e-9; L = np.eye(len(S)) - S / np.sqrt(np.outer(dg, dg)); ev, U = np.linalg.eigh(L)
k = int('${NFAM:-0}') or int(np.argmax(np.diff(ev[:7])[1:]) + 2); E = U[:, :k] / np.maximum(np.linalg.norm(U[:, :k], axis=1)[:, None], 1e-12)
rng = np.random.default_rng(0); best = None
for rep in range(20):   # k-means on the spectral embedding, best of 20 starts
    C = E[rng.choice(len(E), k, replace=False)]
    for it in range(100): lab = np.argmin(((E[:, None] - C[None]) ** 2).sum(-1), 1); C = np.array([E[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
    sse = ((E - C[lab]) ** 2).sum(); best = (sse, lab) if best is None or sse < best[0] else best
lab = best[1]; L_ = np.loadtxt('$R/whole/board/klem_latent.txt', comments='#'); pind, cell = L_[:, 0].astype(int), L_[:, 3].astype(int)
ncell = len(fs); on = (cell >= 1) & (cell <= ncell)
for j in range(k): np.savetxt('$R/expand/family%d_subset.txt' % (j + 1), np.sort(pind[on & (lab[np.clip(cell - 1, 0, ncell - 1)] == j)]), fmt='%d')
open('$R/expand/nfam.txt', 'w').write('%d\n' % k); print('eigenvalues %s -> %d families; %d of %d particles on the grid' % (' '.join('%.3f' % v for v in ev[:7]), k, on.sum(), len(on)))
PYEOF
fi
NF=$(cat $R/expand/nfam.txt); log "$NF families: $(for f in $R/expand/family*_subset.txt; do printf '%s ' $(grep -c . $f); done)"
[ -s $R/periphery.mrc ] || $PY $T/periphery_mask.py $T/maps160/toy_consensus.mrc $R/periphery.mrc --core ${CORE:-0.40} 2>&1 | grep -v Warning | tail -1
[ -s $R/cand/region_1.mrc ] || $PY $T/derive_regions2.py $R/whole/board $R/cand --k ${KREG:-6} --maxreg ${KREG:-6} --within $R/periphery.mrc 2>&1 | grep -v Warning | tail -1   # regions are spatial: derived once
# ---- (2) per family: one focused board per region (the regions of the whole board), 4 at a time x 16 threads
shopt -s nullglob
for F in $(seq 1 $NF); do D=$R/fam$F; S=$R/expand/family${F}_subset.txt; [ -s $S ] || continue; mkdir -p $D/refine/1_refine3D_states; ln -sf $PS $D/refine/1_refine3D_states/p.simple; ln -sf $META $D/meta_fb.npz
  for f in $R/cand/region_*.mrc; do r=$(basename $f .mrc); r=${r#region_}; mkdir -p $D/region_c$r
    while [ $(jobs -rp | wc -l) -ge 4 ]; do sleep 5; done; fb $D/region_c$r/board $f $S "12 10 10 8 8" "1.2 1.0 1.0 0.8 0.8" 16 &
  done; wait; log "family $F ($(grep -c . $S) particles): $(ls -d $D/region_c*/board/klem_latent.txt 2>/dev/null | wc -l) region boards"
done
# ---- (3) readouts: top split vs truth, factors per family with the truth along each motion
$PY $T/org_readout.py $R $T/sim$TAG/gt_particles.txt --fams ${FAMS:-swing:200:1.0} 2>&1 | grep -v Warning | sed -n '/== (1)/,/== (2-4)/p' | head -12
$PY $T/factor_sample.py $R $R/factors --fams $(seq -s, 1 $NF) --glob 'region_c*' --nbin ${NBIN:-8} --budget ${BUDGET:-1000} --minend ${MINEND:-300} --gt $T/sim$TAG/gt_particles.txt --codes ${CODES:?} 2>&1 | grep -v Warning
log DONE

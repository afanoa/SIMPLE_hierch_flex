#!/bin/bash
#SBATCH -J rnasep_focus_vs_focal
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Did the 2026-10-14 FOCUSED FLEX runs (flex_pca pcg_mskfile = region{1,2,3}_focus.mrc, gmm placement, 2 states each,
# ~5-6 % minority on 115k particles) already see the docked state? Each region mask's overlap with the P12/P10.1 focal
# mask, and each run's 2-state labels crossed with the per-particle focal score of the family-1 board (mean score per
# label, docked-top-30k enrichment), after checking that the rows map (defocus / pose agreement).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; FD=${FD:-$A/rnasep_klem/expand1/fam1_pipeline/klem_fb/closed_focal}; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
D=$FD/focus_regions; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') fr] $*"; }
for k in 1 2 3; do f=$A/rnasep_coupled/focus/region$k/1_flex_pca/flex_pca_deconv_labels.txt; echo "-- region$k labels: $(head -3 $f | tr '\n' '|' | cut -c1-200)  ($(grep -c . $f) lines)"; ls $A/rnasep_coupled/focus/region$k/1_flex_pca/*.txt | xargs -n1 basename | tr '\n' ' ' | cut -c1-300; echo; done
ls $A/rnasep_regions/region*_focus.mrc 2>/dev/null || ls $A/rnasep_regions/ | grep -i "focus\|region" | head
$PY - <<PYEOF
import numpy as np, glob, sys, os
sys.path.insert(0, '$W'); from model_class import read_mrc
A = '$A'; M = np.load('$META'); keys = list(M.keys()); print('meta keys:', keys[:30]); N = len(M['flags'])
foc = read_mrc('$FD/focal_mask_p12.mrc'); fm = foc > 0.5
for k in (1, 2, 3):
    fs = glob.glob(A + '/rnasep_regions/region%d_focus.mrc' % k) or glob.glob(A + '/rnasep_regions/*region%d*.mrc' % k)
    if fs:
        rm = read_mrc(fs[0]) > 0.5; print('region %d mask (%s): %d voxels; covers %.0f %% of the P12/P10.1 focal mask, focal mask covers %.0f %% of it' % (k, os.path.basename(fs[0]), rm.sum(), 100. * (rm & fm).sum() / max(fm.sum(), 1), 100. * (rm & fm).sum() / max(rm.sum(), 1)))
F = np.loadtxt('$FD/board/focal_scores.txt', comments='#'); fsc = F[:, 1]; R = np.load('$FD/board/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(R['i1']); s = r[:, :ngrid] @ fsc[:ngrid]
score = np.full(N + 1, np.nan); score[pinds] = s; top = np.zeros(N + 1, bool); top[pinds[np.argsort(-s)[:30000]]] = True
for k in (1, 2, 3):
    f = A + '/rnasep_coupled/focus/region%d/1_flex_pca/flex_pca_deconv_labels.txt' % k
    rows = [l.split() for l in open(f) if l.strip() and not l.startswith('#')]; arr = np.array([[float(x) for x in rr] for rr in rows if len(rr) >= 2])
    pind = arr[:, 0].astype(int); lab = arr[:, 1].astype(int) if arr.shape[1] > 1 else None
    if lab is None or pind.max() > N: print('region %d: unexpected label table (cols %d, max pind %d, N %d)' % (k, arr.shape[1], pind.max(), N)); continue
    ok = ~np.isnan(score[pind]); print('region %d: %d labelled particles, %d on the family-1 board; labels %s' % (k, len(pind), ok.sum(), np.bincount(lab)))
    for l in np.unique(lab):
        m = ok & (lab == l); print('   label %d: %d on the board, mean focal score %+.3f, docked-top-30k fraction %.3f (board average %.3f)' % (l, m.sum(), np.nanmean(score[pind[m]]) if m.any() else float('nan'), top[pind[m]].mean() if m.any() else float('nan'), top[pinds].mean()))
PYEOF
log DONE

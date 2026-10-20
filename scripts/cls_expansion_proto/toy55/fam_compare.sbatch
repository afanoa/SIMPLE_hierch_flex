#!/bin/bash
#SBATCH -J f692_fam_compare
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# family 1 vs family 2 of the real run: view distribution (e1, e2 histograms, overlap), defocus, junk mass on the parent
# board, per-stack share. Readout only.
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan
$PY - <<PY
import numpy as np
M = np.load('$A/rnasep_coupled/basis_span_all/meta_all.npz'); R = np.load('$A/rnasep_klem/f692/parent/klem_resp.npz'); r = R['r']; pinds = R['pinds'].astype(int); ngrid = 49
fam = {f: np.loadtxt('$A/rnasep_klem/f692/expand/family%d_subset.txt' % f, dtype=int) for f in (1, 2)}; pos = {p: i for i, p in enumerate(pinds)}
for f in (1, 2):
    idx = fam[f] - 1; rows = np.array([pos[p] for p in fam[f] if p in pos])
    e1, e2, df = M['e1'][idx], M['e2'][idx], 0.5 * (M['dfx'][idx] + M['dfy'][idx]); junk = r[rows, ngrid:].sum(1) / np.maximum(r[rows].sum(1), 1e-30)
    h2, _, _ = np.histogram2d(e1, e2, bins=[18, 36], range=[[0, 180], [0, 360]]); p = h2 / h2.sum()
    print('family %d: %d particles | e1 quartiles %s | e2 quartiles %s | defocus median %.2f um (q 25-75 %.2f-%.2f) | parent-board junk mass mean %.3f, >0.5: %.3f | top view bins share %.3f' % (f, len(idx), np.percentile(e1, [25, 50, 75]).round(0), np.percentile(e2, [25, 50, 75]).round(0), np.median(df), *np.percentile(df, [25, 75]), junk.mean(), (junk > 0.5).mean(), np.sort(p.ravel())[-10:].sum()))
    if f == 1: p1 = p
    else: print('view-distribution overlap between the families (sum of min of the normalised 2D histograms): %.3f' % np.minimum(p1, p).sum())
stk = M['stkind']
for f in (1, 2):
    s = stk[fam[f] - 1]; print('family %d share per stack: %s' % (f, np.array2string(np.bincount(s, minlength=11)[1:] / len(s), precision=3)))
PY

#!/bin/bash
#SBATCH -J f692_fam_defocus
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# does defocus leak into the latent? parent board: Spearman(z1,z2 ; defocus), mean defocus per argmax cell and per
# family; family boards: mean defocus per k-means state (fam*/states.txt) and per feedback-board cell. Readout only.
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan
$PY - <<PY
import numpy as np
from scipy.stats import spearmanr
M = np.load('$A/rnasep_coupled/basis_span_all/meta_all.npz'); df_all = 0.5 * (M['dfx'] + M['dfy'])
L = np.loadtxt('$A/rnasep_klem/f692/parent/klem_latent.txt', comments='#'); p = L[:, 0].astype(int); z1, z2, am = L[:, 1], L[:, 2], L[:, 3].astype(int); df = df_all[p - 1]
print('parent board: Spearman(z1, defocus) %+.3f  Spearman(z2, defocus) %+.3f' % (spearmanr(z1, df)[0], spearmanr(z2, df)[0]))
cm = np.array([df[am == c].mean() if (am == c).any() else np.nan for c in range(1, 52)]); print('mean defocus per argmax cell (1..49, then free 50,51): %s' % np.array2string(cm, precision=2, max_line_width=250))
fam = {f: np.loadtxt('$A/rnasep_klem/f692/expand/family%d_subset.txt' % f, dtype=int) for f in (1, 2)}
for f in (1, 2):
    d = df_all[fam[f] - 1]; print('family %d: defocus mean %.2f, deciles %s' % (f, d.mean(), np.array2string(np.percentile(d, [10, 30, 50, 70, 90]), precision=2)))
    st = np.loadtxt('$A/rnasep_klem/f692/fam%d/states.txt' % f, dtype=int)
    if len(st) == len(df_all): lab = st[fam[f] - 1]
    else: lab = st
    print('  k-means states of family %d: mean defocus per state %s | state sizes %s' % (f, np.array2string(np.array([df_all[fam[f] - 1][lab == s].mean() for s in range(1, 11)]), precision=2), np.bincount(lab)[1:]))
    Lf = np.loadtxt('$A/rnasep_klem/f692/fam%d/klem_fb/klem_latent.txt' % f, comments='#'); pf = Lf[:, 0].astype(int); dff = df_all[pf - 1]
    print('  feedback board of family %d: Spearman(z1, defocus) %+.3f (z2 %+.3f); cell mean-defocus range %.2f-%.2f' % (f, spearmanr(Lf[:, 1], dff)[0], spearmanr(Lf[:, 2], dff)[0], *np.percentile([dff[Lf[:, 3].astype(int) == c].mean() for c in range(1, 50) if (Lf[:, 3].astype(int) == c).sum() > 100], [0, 100])))
PY

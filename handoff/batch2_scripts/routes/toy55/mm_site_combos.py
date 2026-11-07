#!/usr/bin/env python3
"""Product of the site layer's pooled sites, for a one-pass re-assignment (mm_fassign.sbatch) without any new basis: for each chosen
site k (level r, posterior >= 0.5 = pole) the pole and complement mean maps are the leaf maps averaged by how many of each leaf's
particles fall on that side (even + odd unfiltered halves), Delta_k = pole - complement; combination s (one bit per site) =
mu - sum_k p_k Delta_k + sum_k s_k Delta_k (p_k the pole share, so the population mean is mu). Initial state = the particle's own pole
bits. Writes OUT/comb_NN.mrc, OUT/comb_table.txt (bit 2 = pole), OUT/comb_states.txt.
usage: mm_site_combos.py TREE_DIR OUT SITE[,SITE..] [SMPD of the leaf maps, default 160 x 1.464 / box]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
T, out, sites = sys.argv[1], sys.argv[2], sys.argv[3].split(','); os.makedirs(out, exist_ok=True)
L = [l.split() for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
def loc(p): return p if os.path.exists(p) else os.path.join(T, p[p.index('1_hetero_tree') + len('1_hetero_tree') + 1:])
Mv = np.array([0.5 * (read_mrc(loc(l[5])).astype(np.float64) + read_mrc(loc(l[6])).astype(np.float64)) for l in L]); B = Mv.shape[1]
lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1); N = len(lab); ok = lab > 0
pole = np.array([np.loadtxt('%s/factorial/post/%s.txt' % (T, s), ndmin=2)[:, 0] >= 0.5 for s in sites])   # project rows x site
cnt = np.bincount(lab[ok], minlength=len(L) + 1)[1:].astype(float); mu = (cnt[:, None, None, None] * Mv).sum(0) / cnt.sum()
D, p = [], []
for k, s in enumerate(sites):
    a = np.bincount(lab[ok & pole[k]], minlength=len(L) + 1)[1:].astype(float); b = cnt - a
    P = (a[:, None, None, None] * Mv).sum(0) / a.sum(); C = (b[:, None, None, None] * Mv).sum(0) / b.sum(); D.append(P - C); p.append(a.sum() / cnt.sum())
    print('site %s: pole %d rows (%.3f), |Delta| %.4g' % (s, a.sum(), p[-1], np.linalg.norm(D[-1])))
base = mu - sum(pk * Dk for pk, Dk in zip(p, D)); r = len(sites)
with open(out + '/comb_table.txt', 'w') as ft:
    ft.write('# state ' + ' '.join(sites) + '  (2 = pole)\n')
    for st in range(2 ** r):
        bits = [(st >> k) & 1 for k in range(r)]
        write_mrc('%s/comb_%02d.mrc' % (out, st + 1), (base + sum(bits[k] * D[k] for k in range(r))).astype(np.float32), float(sys.argv[4]) if len(sys.argv) > 4 else 160 * 1.464 / B)
        ft.write('%d %s\n' % (st + 1, ' '.join(str(b_ + 1) for b_ in bits)))
st0 = np.where(ok, 1 + sum(pole[k].astype(int) << k for k in range(r)), 0); rs = np.random.default_rng(5)
for s_ in range(1, 2 ** r + 1):   # refine3D's multi-state startup needs every state populated: 50 random rows of the largest state
    if (st0 == s_).sum() < 50: big = np.bincount(st0[st0 > 0]).argmax(); st0[rs.choice(np.where(st0 == big)[0], 50, replace=False)] = s_
np.savetxt(out + '/comb_states.txt', st0, fmt='%d')
print('%d combination maps (box %d); initial states %s' % (2 ** r, B, ' '.join('%d:%d' % (s_, (st0 == s_).sum()) for s_ in range(1, 2 ** r + 1))))

#!/usr/bin/env python3
"""Truth-expected difference between groups (synthetic sets): every group's expected map = sum over conformers of its truth
composition x the true conformer map (junk: vol_blob / vol_broken by the sim's oris files); per pair the relative difference
energy |D|^2 / |node mean|^2 per band (bands > 40 / 40-20 / 20-10 A) and the composition (class shares, top conformers).
usage: group_truth_diff.py SIMDIR LABELS[,LABELS2] [--maps SIMDIR/maps]"""
import argparse, glob, os, re, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('sim'); ap.add_argument('labels'); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
lab = sum(np.loadtxt(f, dtype=int) for f in a.labels.split(',')); G = np.loadtxt(a.sim + '/gt_particles.txt', comments='#'); conf = G[:, 1].astype(int)[:len(lab)]
jt = np.zeros(len(lab), int)   # junk type: 1 blob, 2 broken (from the oris files' particle counts is not possible per row) -> blob/broken split unknown: use the mean junk map
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
mp = {}
for f in glob.glob(a.sim + '/maps/conf*_m160.mrc'):
    m = re.match(r'conf(\d+)_', os.path.basename(f)); mp[int(m.group(1))] = f
J = 0.5 * (rd(a.sim + '/vol_blob.mrc') + rd(a.sim + '/vol_broken.mrc')) if os.path.isfile(a.sim + '/vol_blob.mrc') else None
cs = sorted(set(conf[lab > 0])); n = rd(mp[next(c for c in cs if c > 0)]).shape[0]
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); res = a.smpd / np.maximum(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2), 1e-9)
bands = [(res > 40), (res <= 40) & (res > 20), (res <= 20) & (res >= 10)]
F = {}
for c in cs: F[c] = np.fft.fftn(rd(mp[c]) if c > 0 else J)
gs = sorted(set(lab[lab > 0])); comp = {g: {c: ((lab == g) & (conf == c)).sum() / max((lab == g).sum(), 1) for c in cs} for g in gs}
M = {g: sum(comp[g][c] * F[c] for c in cs if comp[g][c] > 0) for g in gs}; tot = (lab > 0).sum(); X = sum(M[g] * (lab == g).sum() / tot for g in gs)
def cls(c): return 'junk' if c == 0 else 'swing' if c >= 200 else 'holo' if c >= 101 else 'apo'
for g in gs:
    cc = comp[g]; sh = {k: sum(v for c, v in cc.items() if cls(c) == k) for k in ('apo', 'holo', 'swing', 'junk')}; top = sorted(cc, key=lambda c: -cc[c])[:4]
    print('group %d (n %d): %s | %s' % (g, (lab == g).sum(), ' '.join('%s %.2f' % kv for kv in sh.items()), ' '.join('%d:%.3f' % (c, cc[c]) for c in top)))
for i, p in enumerate(gs):
    for q in gs[i + 1:]:
        D = M[p] - M[q]; print('pair %d-%d: truth relative difference energy per band %s' % (p, q, ' '.join('%.2e' % ((np.abs(D[b]) ** 2).sum() / (np.abs(X[b]) ** 2).sum()) for b in bands)))

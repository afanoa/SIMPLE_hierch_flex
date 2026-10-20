#!/usr/bin/env python3
"""Merge data-derived regions that describe the same moving part: for every pair of region boards (all on the same
particles), the top canonical correlation between their coordinates (quadratic features of z1, z2; fitted on one half of
the particles, measured on the other); pairs above THR are joined (connected components) and each group's masks are
united (voxelwise max). Writes OUT/region_<k>.mrc (largest group first) and OUT/groups.txt. No truth.
--nuis TABLE (ptcl3D print_project_field rows = project rows): every board's features are first regressed on per-particle
nuisances (score, defocus, angular step, shifts; linear + squares), so a shared quality / defocus sorting does not
correlate unrelated regions. THR = auto: median + 3 x 1.4826 x MAD of the off-diagonal correlations (most pairs are
different parts; pairs on one part are the outliers above the bulk). THR = gap: the midpoint of the largest jump between
sorted off-diagonal correlations above their median, used only if the jump exceeds GUARD robust sd of the values below it
(otherwise nothing merges); THR = gap:GUARD sets the guard (default 2).
usage: region_merge.py OUT THR [--nuis TABLE] MASK1:BOARD1 MASK2:BOARD2 ...   (THR a number, auto, gap or gap:GUARD)"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
out, thr = sys.argv[1], sys.argv[2]; args = sys.argv[3:]; nuis = ''
if args and args[0] == '--nuis': nuis, args = args[1], args[2:]
pairs = [a.split(':') for a in args]; os.makedirs(out, exist_ok=True); n = len(pairs)
Z = []
for m, b in pairs:
    L = np.loadtxt(b + '/klem_latent.txt', comments='#'); o = np.argsort(L[:, 0]); Z.append((L[o, 0].astype(int), L[o, 1:3]))
common = Z[0][0]
for p, _ in Z[1:]: common = np.intersect1d(common, p)
def feats(p, z):
    z = z[np.searchsorted(p, common)]; z = (z - z.mean(0)) / z.std(0); return np.c_[z, z ** 2, z[:, :1] * z[:, 1:2]]
F = [feats(p, z) for p, z in Z]
if nuis:
    keys = ('corr', 'dfx', 'dist', 'dist_inpl', 'x', 'y'); T = np.array([[float(dict(w.split('=', 1) for w in l.split() if '=' in w).get(k, 0)) for k in keys] for l in open(nuis) if l.strip() and not l.startswith('#')])
    N_ = T[common - 1]; N_ = (N_ - N_.mean(0)) / np.maximum(N_.std(0), 1e-9); X = np.c_[np.ones(len(N_)), N_, N_ ** 2]
    F = [f - X @ np.linalg.lstsq(X, f, rcond=None)[0] for f in F]; print('nuisances regressed out of every board: %s (linear + squares)' % ' '.join(keys))
half = (np.arange(len(common)) % 2) == 0
def cca1(X, Y):   # top canonical correlation: directions from the even particles, correlation on the odd ones
    Xa, Ya = X[half] - X[half].mean(0), Y[half] - Y[half].mean(0); Cxx, Cyy, Cxy = Xa.T @ Xa, Ya.T @ Ya, Xa.T @ Ya
    Wx = np.linalg.inv(np.linalg.cholesky(Cxx + 1e-9 * np.eye(len(Cxx)))); Wy = np.linalg.inv(np.linalg.cholesky(Cyy + 1e-9 * np.eye(len(Cyy))))
    U, s, Vt = np.linalg.svd(Wx @ Cxy @ Wy.T); a, b = Wx.T @ U[:, 0], Wy.T @ Vt[0]
    return abs(np.corrcoef(X[~half] @ a, Y[~half] @ b)[0, 1])
C = np.eye(n)
for i in range(n):
    for j in range(i + 1, n): C[i, j] = C[j, i] = cca1(F[i], F[j])
print('canonical correlation between region boards (%d common particles; rows = candidates 1..%d):' % (len(common), n))
for i in range(n): print('  %2d %s  %s' % (i + 1, ' '.join('%.2f' % v for v in C[i]), os.path.basename(pairs[i][0])))
off = C[np.triu_indices(n, 1)]
if thr == 'auto': med = np.median(off); thr = med + 3 * 1.4826 * np.median(np.abs(off - med)); print('auto threshold: median %.2f + 3 robust sd -> %.2f' % (med, thr))
if thr.startswith('gap'):
    guard = float(thr.split(':')[1]) if ':' in thr else 2.0
    s = np.sort(off); med = np.median(s); up = np.where(s >= med)[0]; g = np.diff(s[up]); k = int(np.argmax(g)) if len(g) else 0
    lo = s[:up[k] + 1] if len(g) else s; sd = 1.4826 * np.median(np.abs(lo - np.median(lo)))
    thr = 0.5 * (s[up[k]] + s[up[k] + 1]) if len(g) and g[k] > guard * sd else 2.0
    print('gap threshold: largest jump %.2f -> %.2f (robust sd below %.3f, guard %.1f sd) -> %s' % (s[up[k]] if len(g) else 0, s[up[k] + 1] if len(g) else 0, sd, guard, '%.2f' % thr if thr < 2 else 'no merge'))
thr = float(thr); lab = list(range(n))
def find(i):
    while lab[i] != i: i = lab[i]
    return i
for i in range(n):
    for j in range(i + 1, n):
        if C[i, j] > thr: lab[find(j)] = find(i)
groups = {}
for i in range(n): groups.setdefault(find(i), []).append(i)
masks = [read_mrc(m) for m, _ in pairs]; G = sorted(groups.values(), key=lambda g: -sum(masks[i].sum() for i in g))
with open(out + '/groups.txt', 'w') as f:
    for k, g in enumerate(G, 1):
        M = np.max([masks[i] for i in g], 0); write_mrc('%s/region_%d.mrc' % (out, k), M.astype(np.float32), 1.464)
        f.write('%d %s\n' % (k, ' '.join(pairs[i][0] for i in g))); print('merged region %d <- candidates %s' % (k, ' '.join(str(i + 1) for i in g)))
print('%d candidates -> %d merged regions (threshold %.2f)' % (n, len(G), thr))

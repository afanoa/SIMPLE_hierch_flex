#!/usr/bin/env python3
"""Is a node's conformational landscape reproducible across independent runs even when its partitions are not? Per run: its K leaf maps
(stage-B even/odd unfiltered half-maps), low-passed (LP A), weighted by the node mask (rsmask.mrc, Fourier-cropped to the map box);
deviation maps d_k = m_k - sum_k p_k m_k (p_k = leaf fractions). Per map: half-map reliability R_k (Spearman-Brown of the even/odd
deviation correlation). Across runs A, B: (1) each map of A explained by B's maps (least-squares fit of d_k^A on all d_j^B: R^2 / R_k^A,
1 = the map lies in B's landscape as far as its noise allows), best single-map match corr / sqrt(R R); (2) principal angles between the
runs' top-3 leaf-size-weighted deviation subspaces; (3) occupancy: leaf-size-weighted distribution of each run's maps on the pooled top-2
PCs (mean, sd, 5-95 % range) and the 1-D Wasserstein distance along PC1.  usage: landscape_cmp.py NAME=MAPDIR:LABELS ... [--lp 15]"""
import sys, itertools, numpy as np
args = [a for a in sys.argv[1:] if '=' in a and not a.startswith('--')]; LP = float(sys.argv[sys.argv.index('--lp') + 1]) if '--lp' in sys.argv else 15.
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def cell(p): return float(np.fromfile(p, np.float32, 13)[10])
def crop(v, n):   # Fourier crop/pad a box-N map to box n at the same physical size
    N = v.shape[0]
    if N == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); lo = N // 2 - n // 2; G = F[lo:lo + n, lo:lo + n, lo:lo + n]
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / N) ** 3
def lowpass(v, smpd, lp):
    n = v.shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij')
    return np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (lp / 1.2) ** 2)))
RUNS = {}
for a in args:
    name, rest = a.split('=', 1); mdir, labf = rest.split(':'); lab = np.loadtxt(labf, dtype=int); ks = sorted(set(lab[lab > 0]))
    E, O = [], []
    for k in ks:
        b = '%s/recvol_state%02d' % (mdir, k); E.append(rd(b + '_even_unfil.mrc')); O.append(rd(b + '_odd_unfil.mrc'))
    n = E[0].shape[0]; smpd = cell(b + '_even_unfil.mrc') / n; RUNS[name] = dict(E=E, O=O, n=n, smpd=smpd, p=np.array([(lab == k).sum() for k in ks], float), ks=ks)
n = RUNS[list(RUNS)[0]]['n']; smpd = RUNS[list(RUNS)[0]]['smpd']
M = crop(rd(K_ + '/rflex_f692st/node_r1/rsmask.mrc'), n); M = np.clip(M, 0, 1); w = M.ravel() / M.sum()
def prep(R):
    p = R['p'] / R['p'].sum(); E = np.array([lowpass(x, R['smpd'], LP).ravel() for x in R['E']]); O = np.array([lowpass(x, R['smpd'], LP).ravel() for x in R['O']])
    F = 0.5 * (E + O); dF, dE, dO = (X - (p[:, None] * X).sum(0) for X in (F, E, O))
    wc = lambda a, b: (w * a * b).sum() / np.sqrt((w * a * a).sum() * (w * b * b).sum())
    r = np.array([wc(dE[k], dO[k]) for k in range(len(p))]); R2 = 2 * r / (1 + r)
    return dict(p=p, d=dF * np.sqrt(w)[None], rel=R2)
P = {k: prep(v) for k, v in RUNS.items()}
for k, v in P.items():
    print('%-7s %d leaves; leaf fractions %s; half-map reliability of the deviation maps %s' % (k, len(v['p']), ' '.join('%.2f' % x for x in v['p']), ' '.join('%.2f' % x for x in v['rel'])))
print('== (1) each map of run A explained by run B\'s maps: fit R^2 / own reliability (1 = within noise), and best single match / noise ceiling')
for a, b in itertools.permutations(P, 2):
    A_, B_ = P[a], P[b]; X = B_['d'].T; out = []
    for k in range(len(A_['p'])):
        y = A_['d'][k]; c, *_ = np.linalg.lstsq(X, y, rcond=None); r2 = 1 - ((y - X @ c) ** 2).sum() / (y ** 2).sum()
        cc = [(y * B_['d'][j]).sum() / np.sqrt((y * y).sum() * (B_['d'][j] ** 2).sum()) / np.sqrt(max(A_['rel'][k] * B_['rel'][j], 1e-6)) for j in range(len(B_['p']))]
        out.append('%.2f/%.2f' % (r2 / max(A_['rel'][k], 1e-6), max(cc)))
    print('  %s by %s: %s' % (a, b, '  '.join(out)))
def subspace(v, m=3):
    D = v['d'] * np.sqrt(v['p'])[:, None]; U, s, Vt = np.linalg.svd(D, full_matrices=False); return Vt[:m], (s ** 2 / (s ** 2).sum())[:m]
print('== (2) principal angles (cosines) between the runs\' top-3 weighted deviation subspaces')
for a, b in itertools.combinations(P, 2):
    Va, ea = subspace(P[a]); Vb, eb = subspace(P[b]); cs = np.linalg.svd(Va @ Vb.T, compute_uv=False)
    print('  %s vs %s: %s   (variance shares %s | %s)' % (a, b, ' '.join('%.2f' % x for x in cs), ' '.join('%.2f' % x for x in ea), ' '.join('%.2f' % x for x in eb)))
print('== (3) occupancy on the pooled top-2 PCs (leaf-size weighted): mean / sd / 5-95 %% range, and the PC1 Wasserstein distance')
allD = np.vstack([v['d'] * np.sqrt(v['p'])[:, None] for v in P.values()]); U, s, Vt = np.linalg.svd(allD, full_matrices=False); B2 = Vt[:2]
def wq(x, p, q):
    o = np.argsort(x); c = np.cumsum(p[o]); return np.interp(q, c - p[o] / 2, x[o])
coords = {}
for k, v in P.items():
    C = v['d'] @ B2.T; coords[k] = C; s_ = []
    for j in range(2):
        m = (v['p'] * C[:, j]).sum(); sd = np.sqrt((v['p'] * (C[:, j] - m) ** 2).sum()); s_.append('PC%d %.3f / %.3f / %.3f..%.3f' % (j + 1, m, sd, wq(C[:, j], v['p'], 0.05), wq(C[:, j], v['p'], 0.95)))
    print('  %-7s %s' % (k, ' | '.join(s_)))
qs = np.linspace(0.01, 0.99, 99)
for a, b in itertools.combinations(P, 2):
    wa = np.mean(np.abs(wq(coords[a][:, 0], P[a]['p'], qs) - wq(coords[b][:, 0], P[b]['p'], qs))); spread = np.sqrt((P[a]['p'] * coords[a][:, 0] ** 2).sum())
    print('  %s vs %s: PC1 Wasserstein %.3f (= %.2f of run %s\'s PC1 spread)' % (a, b, wa, wa / max(spread, 1e-12), a))

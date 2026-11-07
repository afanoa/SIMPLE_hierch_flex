#!/usr/bin/env python3
"""lsamp_repro_cmp2.py = lsamp_repro_cmp.py with every run after the first oriented to correlate positively with the first run's full
coordinate (the first run by its arm readout), for axes that are not the closed-open direction. lsamp_repro_cmp.py: cross-run reproducibility of one node's landscape (lsamp_repro.sbatch outputs of several independent tree runs on
the same particles). Per run: axes, sample count, closed-arm / open-arm of every sample (held-out maps, both folds; own-pixel crop
maps), occupancy (sample sizes). Across runs: per-particle coordinate correlation (common rows; each run oriented so its first sample
is the closed end), sample-label agreement (NMI, same or adjacent sample), end-sample particle overlap (Jaccard), and matched maps:
held-out maps of run X (pixel set B) against run Y's (pixel set A) - independent noise - correlated at lp 15 inside the sphere.
usage: lsamp_repro_cmp2.py OUTFILE NAME:DIR [NAME:DIR ...]"""
import sys, glob, os, itertools, numpy as np
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def arm(f):
    out = {}
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        w = l.split(); k = int(os.path.basename(w[0])[6:8]); out[k] = (float(w[1]), float(w[2]))
    return out
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; return (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum() / max(np.sqrt((px * np.log(px)).sum() * (py * np.log(py)).sum()), 1e-12)
outf = sys.argv[1]; runs = [a.split(':', 1) for a in sys.argv[2:]]; R = {}; lines = []
def P(*s): line = ' '.join(str(x) for x in s); print(line, flush=True); lines.append(line)
n = 160; smpd = 1.464; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=smpd), indexing='ij'); G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (15. / 1.2) ** 2)
g = np.arange(n) - n // 2; sph = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) < 0.48 * n
def prep(fn): v = np.fft.irfftn(np.fft.rfftn(rd(fn)) * G, s=(n, n, n))[sph]; v -= v.mean(); return v / np.linalg.norm(v)
for name, d in runs:
    lab = np.loadtxt(d + '/samples.txt', dtype=int); C = np.load(d + '/coords.npz'); ah = [arm(d + '/arm_heldout_%s.txt' % f_) for f_ in 'AB']; af = arm(d + '/arm_final_all.txt')
    ns = lab.max(); sizes = np.bincount(lab, minlength=ns + 1)[1:]
    if not R: flip = af[1][0] < af[ns][0]; t0 = np.full(len(lab), np.nan); t0[C['rows']] = C['tfull'][:, 0] * (-1 if flip else 1)   # first run: sample 1 = the closed end
    else: tt = np.full(len(lab), np.nan); tt[C['rows']] = C['tfull'][:, 0]; mm_ = np.isfinite(tt) & np.isfinite(t0); flip = np.corrcoef(tt[mm_], t0[mm_])[0, 1] < 0   # the others follow it
    sg = -1 if flip else 1; t = np.full(len(lab), np.nan); t[C['rows']] = C['tfull'][:, 0] * sg
    tAf = np.full(len(lab), np.nan); tAf[C['rows']] = C['tA'][:, 0] * sg; tBf = np.full(len(lab), np.nan); tBf[C['rows']] = C['tB'][:, 0] * sg
    o = (lambda k: ns + 1 - k) if flip else (lambda k: k); labo = np.where(lab > 0, np.array([0] + [o(k) for k in range(1, ns + 1)])[lab], 0)
    fo = {}   # the fold cells have their own count (sparse merges may differ): orient each fold by the same flip
    for f_ in 'AB':
        fl = sorted(glob.glob(d + '/heldout_%s/state_*.mrc' % f_)); nf = len(fl); fo[f_] = (lambda k, nf=nf: nf + 1 - k) if flip else (lambda k: k)
    R[name] = dict(d=d, lab=labo, t=t, ns=ns, sizes=[sizes[o(k) - 1] for k in range(1, ns + 1)], ah=[{fo[f_](k): v for k, v in a.items()} for f_, a in zip('AB', ah)], af={o(k): v for k, v in af.items()},
                   heldA={fo['A'](int(os.path.basename(p)[6:8])): p for p in glob.glob(d + '/heldout_A/state_*.mrc')}, heldB={fo['B'](int(os.path.basename(p)[6:8])): p for p in glob.glob(d + '/heldout_B/state_*.mrc')},
                   nper=C['nper'], ax=C['ax'], tA=tAf, tB=tBf)
    r = R[name]; mm = np.isfinite(r['tA']) & np.isfinite(r['tB']); r['rAB'] = np.corrcoef(r['tA'][mm], r['tB'][mm])[0, 1]
    P('== %s: axes %s, %s samples (closed end first); sizes %s; fold coordinates A vs B corr %.3f (independent pixel noise: the signal part)' % (name, ' '.join(str(j + 1) for j in r['ax']), 'x'.join(map(str, r['nper'])), ' '.join(map(str, r['sizes'])), r['rAB']))
    nk = max(len(r['heldA']), len(r['heldB']))
    P('   held-out closed-arm (fold A / B): %s' % '  '.join('%d:%.2f/%.2f' % (k, r['ah'][0].get(k, (np.nan,))[0], r['ah'][1].get(k, (np.nan,))[0]) for k in range(1, nk + 1)))
    P('   held-out open-arm   (fold A / B): %s' % '  '.join('%d:%.2f/%.2f' % (k, r['ah'][0].get(k, (np.nan, np.nan))[1], r['ah'][1].get(k, (np.nan, np.nan))[1]) for k in range(1, nk + 1)))
    P('   own-pixel closed-arm (final samples, crop maps): %s' % '  '.join('%d:%.2f' % (k, r['af'][k][0]) for k in range(1, ns + 1)))
P('== across runs (same %d particles)' % int(np.isfinite(R[runs[0][0]]['t']).sum()))
for (a, ra), (b, rb) in itertools.combinations(R.items(), 2):
    m = np.isfinite(ra['t']) & np.isfinite(rb['t']) & (ra['lab'] > 0) & (rb['lab'] > 0)
    rc = np.corrcoef(ra['t'][m], rb['t'][m])[0, 1]
    la, lb = ra['lab'][m], rb['lab'][m]; fa, fb = (la - 1) / max(ra['ns'] - 1, 1), (lb - 1) / max(rb['ns'] - 1, 1)   # sample positions on [0, 1]
    near = np.mean(np.abs(fa - fb) <= 1.0 / max(min(ra['ns'], rb['ns']) - 1, 1) + 1e-9)
    ja = lambda k1, k2: np.sum((la == k1) & (lb == k2)) / max(np.sum((la == k1) | (lb == k2)), 1)
    mx = np.isfinite(ra['tA']) & np.isfinite(rb['tB']) & np.isfinite(ra['tB']) & np.isfinite(rb['tA'])
    cr = 0.5 * (np.corrcoef(ra['tA'][mx], rb['tB'][mx])[0, 1] + np.corrcoef(ra['tB'][mx], rb['tA'][mx])[0, 1]); ag = cr / np.sqrt(max(ra['rAB'] * rb['rAB'], 1e-12))
    P('%s vs %s: coordinate corr %.3f (cross-pixel %.3f = %.2f of the within-run ceiling); sample NMI %.3f; same-or-adjacent position %.1f %%; closed-end Jaccard %.2f, open-end Jaccard %.2f' % (a, b, rc, cr, ag, nmi(la, lb), 100 * near, ja(1, 1), ja(ra['ns'], rb['ns'])))
    # matched held-out maps: run a fold A (pixels B) vs run b fold B (pixels A), every pair of samples; best match per sample of a
    ka, kb = sorted(ra['heldA']), sorted(rb['heldB']); VA = {k: prep(ra['heldA'][k]) for k in ka}; VB = {k: prep(rb['heldB'][k]) for k in kb}
    Mx = np.array([[float(VA[i] @ VB[j]) for j in kb] for i in ka]); best = [kb[int(np.argmax(row))] for row in Mx]
    P('   matched held-out maps (%s fold A vs %s fold B, lp 15, corr): per %s sample the best %s sample: %s' % (a, b, a, b, ' '.join('%d->%d(%.3f)' % (i, j, Mx[ii].max()) for ii, (i, j) in enumerate(zip(ka, best)))))
    P('   ends: %s closed end vs %s closed end corr %.3f, open ends %.3f; closed vs open within %s %.3f' % (a, b, Mx[0, 0], Mx[-1, -1], a, Mx[0, -1]))
open(outf, 'w').write('\n'.join(lines) + '\n')

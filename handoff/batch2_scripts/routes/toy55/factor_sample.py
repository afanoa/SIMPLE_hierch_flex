#!/usr/bin/env python3
"""Factor sampling: continuous heterogeneity as factors, one per independently moving part, each sampled evenly end to end.
Per family of RUN: every region board's per-particle coordinates (klem_latent z1 z2, rows = project rows); boards whose
coordinates depend on each other (cross-validated canonical correlation of quadratic features, fitted on one half of the
particles and measured on the other, >= RMTHR) are one factor, represented by its most confident board (mean max responsibility);
the factor's coordinate t = that board's main axis (confidence-weighted PCA of z1 z2; the second axis's variance ratio flags a
sheet); NBIN equal-width bins between the 1st and 99th percentile of t (end bins open to the tails); an interior bin below BUDGET
joins its smaller neighbour, an end bin below MINEND joins inward. Each factor is written as a state set (OUT/fam<F>_factor<N>/
states.txt over project rows + groups.txt). GT + CODES: per bin the mean +- sd of every truth coordinate (codes name:id0:step
for one motion per conformer range, grid:id0:kstep:jstep for conformer id0 + 10 k + j carrying two motions).
usage: factor_sample.py RUN OUT [--fams 1,2 --glob region* --rmthr 0.5 --nbin 10 --budget 3000 --minend 1000 --gt GT --codes ...]"""
import sys, os, glob, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('out'); ap.add_argument('--fams', default='1,2'); ap.add_argument('--glob', default='region*')
ap.add_argument('--rmthr', type=float, default=0.5); ap.add_argument('--nbin', type=int, default=10); ap.add_argument('--budget', type=int, default=3000); ap.add_argument('--minend', type=int, default=1000)
ap.add_argument('--gt', default=''); ap.add_argument('--codes', default=''); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
truth = {}; mf = '%s/fam1/meta_fb.npz' % a.run; nrow = len(np.load(mf)['flags']) if os.path.isfile(mf) else None   # project rows
if a.gt:
    G = np.loadtxt(a.gt, comments='#'); conf = np.zeros(int(G[:, 0].max()) + 1, int); conf[G[:, 0].astype(int)] = G[:, 1].astype(int); nrow = nrow or len(conf) - 1
    for c in [w for w in a.codes.split(',') if w]:
        p = c.split(':')
        if p[0] == 'table':   # table:FILE:name - a per-conformer value (e.g. the closed-arm occupancy of its map) for every particle
            tv = {int(l.split()[0]): float(l.split()[1]) for l in open(p[1]) if l.strip() and not l.startswith('#')}; truth[p[2]] = np.array([tv.get(int(c), np.nan) for c in conf]); continue
        if p[0] == 'grid':
            i0, ks, js = int(p[1]), float(p[2]), float(p[3]); m = conf >= i0; truth['grid_k'] = np.where(m, ((conf - i0) // 10) * ks, np.nan); truth['grid_j'] = np.where(m, ((conf - i0) % 10) * js, np.nan)
        else:
            nm, i0, st = p[0], int(p[1]), float(p[2]); hi = min([int(q.split(':')[1]) for q in a.codes.split(',') if q and not q.startswith(('table:', 'grid:')) and int(q.split(':')[1]) > i0] + [10 ** 9])
            truth[nm] = np.where((conf >= i0) & (conf < hi), (conf - i0) * st, np.nan)
def feats(Z): z1, z2 = Z[:, 0], Z[:, 1]; F = np.c_[z1, z2, z1 ** 2, z2 ** 2, z1 * z2]; return F
def cca_cv(X, Y):   # first canonical correlation, fitted on even rows and measured on odd rows and vice versa
    r = []
    for f in (0, 1):
        tr, te = (np.arange(len(X)) % 2) == f, (np.arange(len(X)) % 2) != f
        mx, my = X[tr].mean(0), Y[tr].mean(0); sx, sy = X[tr].std(0) + 1e-9, Y[tr].std(0) + 1e-9; A, B = (X[tr] - mx) / sx, (Y[tr] - my) / sy
        Cxx, Cyy, Cxy = A.T @ A / len(A) + 1e-6 * np.eye(A.shape[1]), B.T @ B / len(B) + 1e-6 * np.eye(B.shape[1]), A.T @ B / len(A)
        ex, Ux = np.linalg.eigh(Cxx); ey, Uy = np.linalg.eigh(Cyy); Wx = Ux @ np.diag(ex ** -0.5) @ Ux.T; Wy = Uy @ np.diag(ey ** -0.5) @ Uy.T
        U, S, Vt = np.linalg.svd(Wx @ Cxy @ Wy); u, v = Wx @ U[:, 0], Wy @ Vt[0]
        p, q = ((X[te] - mx) / sx) @ u, ((Y[te] - my) / sy) @ v; r.append(abs(np.corrcoef(p, q)[0, 1]))
    return float(np.mean(r))
for F in [int(x) for x in a.fams.split(',')]:
    B = {}
    for f in sorted(glob.glob('%s/fam%d/%s/board/klem_latent.txt' % (a.run, F, a.glob))):
        L = np.loadtxt(f, comments='#'); B[f.split('/')[-3]] = (L[:, 0].astype(int), L[:, 1:3], L[:, 4])
    if not B: print('family %d: no boards' % F); continue
    common = sorted(set.intersection(*[set(v[0]) for v in B.values()])); names = list(B); idx = {}
    for b in names:
        pos = {p: i for i, p in enumerate(B[b][0])}; idx[b] = np.array([pos[p] for p in common])
    rows = np.array(common); Fe = {b: feats(B[b][1][idx[b]]) for b in names}; nb = len(names); Cm = np.eye(nb)
    for i in range(nb):
        for j in range(i + 1, nb): Cm[i, j] = Cm[j, i] = cca_cv(Fe[names[i]], Fe[names[j]])
    print('== family %d: %d region boards on %d common particles; cross-validated canonical correlation between their coordinates:' % (F, nb, len(rows)))
    for i in range(nb): print('   %-16s %s' % (names[i], ' '.join('%.2f' % v for v in Cm[i])))
    lab = list(range(nb))
    for i in range(nb):
        for j in range(nb):
            if Cm[i, j] >= a.rmthr: li, lj = lab[i], lab[j]; lab = [li if l == lj else l for l in lab]
    groups = {}
    for i, l in enumerate(lab): groups.setdefault(l, []).append(names[i])
    groups = sorted(groups.values(), key=lambda g: -len(g)); print('   %d factors at %.2f: %s' % (len(groups), a.rmthr, ' | '.join('+'.join(g) for g in groups)))
    for fi, g in enumerate(groups, 1):
        rep = max(g, key=lambda b: B[b][2].mean()); Z = B[rep][1][idx[rep]]; w = B[rep][2][idx[rep]]; mu = (w[:, None] * Z).sum(0) / w.sum()
        C = ((Z - mu) * w[:, None]).T @ (Z - mu) / w.sum(); ev, V = np.linalg.eigh(C); t = (Z - mu) @ V[:, -1]; ratio = ev[0] / max(ev[-1], 1e-12)
        q1, q99 = np.percentile(t, [1, 99]); edges = np.linspace(q1, q99, a.nbin + 1); bins = np.clip(np.searchsorted(edges, t, side='right') - 1, 0, a.nbin - 1)
        segs = [list(np.where(bins == k)[0]) for k in range(a.nbin)]; segs = [s for s in segs if s]
        i = 1
        while i < len(segs) - 1:
            if len(segs[i]) < a.budget:
                jn = i - 1 if len(segs[i - 1]) <= len(segs[i + 1]) else i + 1; segs[jn] = segs[jn] + segs[i]; segs.pop(i); i = max(1, i - 1)
            else: i += 1
        for e in (0, -1):
            if len(segs) > 2 and len(segs[e]) < a.minend: k = 1 if e == 0 else -2; segs[k] = segs[k] + segs[e]; segs.pop(e)
        d = '%s/fam%d_factor%d' % (a.out, F, fi); os.makedirs(d, exist_ok=True); n_out = nrow or int(rows.max())
        st = np.zeros(n_out, int)
        for k, s in enumerate(segs, 1): st[rows[s] - 1] = k
        np.savetxt(d + '/states.txt', st, fmt='%d')
        with open(d + '/groups.txt', 'w') as fo:
            fo.write('# factor %d of family %d (boards %s, coordinate from %s, second-axis variance ratio %.2f): bin, particles, name\n' % (fi, F, '+'.join(g), rep, ratio))
            for k, s in enumerate(segs, 1): fo.write('%d %d f%d_factor%d_bin%02d\n' % (k, len(s), F, fi, k))
        line = '   factor %d (%s; coordinate %s, sheet ratio %.2f): %d bins, particles %s' % (fi, '+'.join(g), rep, ratio, len(segs), ' '.join(str(len(s)) for s in segs))
        print(line)
        t2 = (Z - mu) @ V[:, 0]
        def rs(x, y):   # Spearman rank correlation where the truth is defined
            m = np.isfinite(y); rx, ry = np.argsort(np.argsort(x[m])), np.argsort(np.argsort(y[m])); return float(np.corrcoef(rx, ry)[0, 1]) if m.sum() > 50 else float('nan')
        if truth: print('      rank correlation with the first / second axis: %s' % '  '.join('%s %+.2f / %+.2f' % (nm, rs(t, tv[rows]), rs(t2, tv[rows])) for nm, tv in truth.items()))
        for nm, tv in truth.items():
            vals = [tv[rows[s]] for s in segs]; mv = [np.nanmean(v) if np.isfinite(v).any() else np.nan for v in vals]
            fr = [np.mean(np.isfinite(v)) for v in vals]
            print('      truth %-8s per bin mean (share of bin): %s' % (nm, ' '.join('%.3g(%.0f%%)' % (m, 100 * s) for m, s in zip(mv, fr))))

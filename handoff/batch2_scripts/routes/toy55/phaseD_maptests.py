#!/usr/bin/env python3
"""usage: phaseD_maptests.py SERIES_DIR:MASK:MSKDIAM[:TRUTHKEY] ..   (phaseC_maptests.py without truth, plus where the change sits)
Without TRUTHKEY there is no truth-implied beta. It also prints the location of the change, from the raw low-passed P5 - P1:
- the focus fraction (the share of the root mask above 3 sd of its half-map noise);
- the centroids (A from the box centre; radius from and height along z, the symmetry axis) of the top 1 % abs(P5 - P1)
  voxels and of the top 0.5 % positive and negative lobes, and the lobe separation.
The original phaseC description follows.

phaseB_maptests.py's step profile and end-map mixture test in three variants, at 15 and 10 A:
- A: raw maps in the root mask;
- B: amplitude-matched maps in the root mask;
- C: amplitude-matched maps in a focus mask of the moving region.
Amplitude matching: per Fourier shell, each bin's signal power is the cross power of its half maps, Re<E O*>. All three of its maps
(full, even, odd) are scaled by sqrt(mean signal power over the bins / its own), clipped to [0.2, 5].
Focus mask: the root mask where the low-passed matched P5 - P1 exceeds 3 x the sd of its half-map noise ((E5 - O5) - (E1 - O1)) / 2.
Also a per-shell mixture fit:
- S: root mask; SF: a focus mask from the RAW low-passed P5 - P1 (same 3 sd rule).
- Each map is multiplied by the soft mask and Fourier transformed. Per shell up to the test resolution, M_k = a(s) P1 + b(s) P5 is fit
  with real coefficients, which absorbs per-bin blur / contrast. The residual power is compared with the half-map noise expectation,
  the cross-fitted half residuals give rho_res, and beta = the power-weighted mean of b / (a + b).
The beta implied by the truth is (t_k - t_1) / (t_5 - t_1), from bins.txt's per-bin open / A / closed fraction or angle mean."""
import sys, re, numpy as np
def read_mrc(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = h.view(np.float32); n = int(h[0])
    return np.frombuffer(raw[1024 + int(h[23]):], np.float32, count=n ** 3).astype(np.float64).reshape(n, n, n), float(hf[10]) / n
def fresample(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2
    if n < m: s = c - n // 2; G = F[s:s + n, s:s + n, s:s + n]
    else: G = np.zeros((n, n, n), complex); s = n // 2 - c; G[s:s + m, s:s + m, s:s + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / m) ** 3
def kr(n): f = np.fft.fftfreq(n) * n; return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def lp_w(n, smpd, res):
    r = kr(n) / n; fc = smpd / res; w = np.clip((1.1 * fc - r) / (0.2 * fc), 0., 1.); return 0.5 - 0.5 * np.cos(np.pi * w)
def corr(a, b): a = a - a.mean(); b = b - b.mean(); return float(a @ b / np.sqrt((a @ a) * (b @ b)))
def truth(D, key):
    t = []
    for l in open(D + '/bins.txt'):
        if not l.startswith('bin '): continue
        m = re.search(r'%s mean ([+-]?[\d.]+)' % 'a\\d_\\w\\w', l) if key == 'mean' else re.search(r'\b%s ([\d.]+)' % key, l); t.append(float(m.group(1)))
    return np.array(t)
for spec in sys.argv[1:]:
    parts = spec.split(':'); D, MK, MD = parts[:3]; KEY = parts[3] if len(parts) > 3 else None; F = {}
    for k in range(1, 6):
        F[k] = [np.fft.fftn(read_mrc(f)[0]) for f in ('%s/bin%d.mrc' % (D, k), '%s/halves/bin%d_even.mrc' % (D, k), '%s/halves/bin%d_odd.mrc' % (D, k))]
    smpd = read_mrc('%s/bin1.mrc' % D)[1]; n = F[1][0].shape[0]; mk = fresample(read_mrc(MK)[0], n) > 0.5; sh = np.rint(kr(n)).astype(int).ravel(); ns = sh.max() + 1
    sp = {k: np.bincount(sh, (F[k][1] * F[k][2].conj()).real.ravel(), ns) / np.maximum(np.bincount(sh, None, ns), 1) for k in F}; spm = np.mean([sp[k] for k in sp], 0)
    scale = {k: np.clip(np.sqrt(np.maximum(spm, 1e-30) / np.maximum(sp[k], 1e-30 + 1e-6 * np.abs(spm))), 0.2, 5.)[sh].reshape(n, n, n) for k in sp}
    t = truth(D, KEY) if KEY else np.full(5, np.nan); bimp = (t[1:4] - t[0]) / (t[4] - t[0]) if KEY else np.full(3, np.nan)
    print('######## %s (box %d, %.3f A/px, root mask %d voxels; truth %s per bin %s; truth-implied beta bins 2-4 %s)' % (D, n, smpd, mk.sum(), KEY, ' '.join('%.3f' % x for x in t), ' '.join('%.2f' % x for x in bimp)))
    w15 = lp_w(n, smpd, 15.); Dr = np.real(np.fft.ifftn((F[5][0] - F[1][0]) * w15)); Ndr = np.real(np.fft.ifftn(((F[5][1] - F[5][2]) - (F[1][1] - F[1][2])) * w15)) / 2
    focf = (mk & (np.abs(Dr) > 3 * Ndr[mk].std())).sum() / mk.sum(); g = np.arange(n) - n // 2; ZZ, YY, XX = np.meshgrid(g, g, g, indexing='ij')
    def cen(sel): z_, y_, x_ = ZZ[sel].mean() * smpd, YY[sel].mean() * smpd, XX[sel].mean() * smpd; return x_, y_, z_, np.hypot(x_, y_)
    ad = np.abs(Dr); top = mk & (ad >= np.percentile(ad[mk], 99)); pos = mk & (Dr >= np.percentile(Dr[mk], 99.5)); neg = mk & (Dr <= np.percentile(Dr[mk], 0.5))
    ct, cp, cn = cen(top), cen(pos), cen(neg); rr = np.sqrt((XX[top] * smpd) ** 2 + (YY[top] * smpd) ** 2)
    print('  where (raw P5 - P1 at 15 A): focus fraction %.3f of the mask; top 1 %% abs centroid x %.0f y %.0f z %.0f A (radius from z %.0f, median voxel radius %.0f A of mask radius %.0f); + lobe (%.0f, %.0f, %.0f) - lobe (%.0f, %.0f, %.0f), separation %.0f A; + / - power ratio %.2f'
          % (focf, ct[0], ct[1], ct[2], ct[3], np.median(rr), float(MD) / 2, cp[0], cp[1], cp[2], cn[0], cn[1], cn[2], np.sqrt(sum((cp[i] - cn[i]) ** 2 for i in range(3))), (Dr[mk] ** 2 * (Dr[mk] > 0)).sum() / max((Dr[mk] ** 2 * (Dr[mk] < 0)).sum(), 1e-30)))
    for res in (15., 10.):
        w = lp_w(n, smpd, res)
        for var in ('A',):
            M = {k: [np.real(np.fft.ifftn(x * w * (scale[k] if var != 'A' else 1.))) for x in F[k]] for k in F}
            if var == 'C':
                Dm = M[5][0] - M[1][0]; Nd = ((M[5][1] - M[5][2]) - (M[1][1] - M[1][2])) / 2; foc = mk & (np.abs(Dm) > 3 * Nd[mk].std())
                if foc.sum() < 0.005 * mk.sum(): foc = mk & (np.abs(Dm) >= np.percentile(np.abs(Dm[mk]), 95))
                use = foc
            else: use = mk
            V = {k: [x[use] for x in M[k]] for k in M}; neff = use.sum() * smpd ** 3 / (res / 2) ** 3
            eo = {k: 1 - corr(V[k][1], V[k][2]) for k in V}; d = [1 - corr(V[k][0], V[k + 1][0]) for k in range(1, 5)]; rn = [d[k - 1] / ((eo[k] + eo[k + 1]) / 4) for k in range(1, 5)]
            P1, P5 = V[1][0], V[5][0]; nv = lambda k: np.var(V[k][1] - V[k][2]) / 4; out = []
            for k in (2, 3, 4):
                A = np.column_stack([P1, P5, np.ones(len(P1))]); a, b, c = np.linalg.lstsq(A, V[k][0], rcond=None)[0]; rr = V[k][0] - A @ np.array([a, b, c])
                exp = nv(k) + a * a * nv(1) + b * b * nv(5); re_ = V[k][1] - (a * V[1][2] + b * V[5][2] + c); ro = V[k][2] - (a * V[1][1] + b * V[5][1] + c); rho = corr(re_, ro)
                beta = np.linalg.lstsq(np.column_stack([P5 - P1, np.ones(len(P1))]), V[k][0] - P1, rcond=None)[0][0]
                out.append('bin%d %.2f / %.3f (z %.1f) / b %.2f' % (k, np.var(rr) / exp, rho, rho * np.sqrt(neff), beta))
            print('  %2.0f A %s%s: steps %s, dominant share %.2f | mixture resid/noise / rho_res (z) / beta: %s' % (res, var, ' (focus %.3f of mask)' % (use.sum() / mk.sum()) if var == 'C' else '',
                  ' '.join('%.1f' % x for x in rn), max(rn) / sum(rn), ' | '.join(out)))
        for var in ('S', 'SF'):
            if var == 'S': use = mk.astype(float)
            else:
                Mr = {k: [np.real(np.fft.ifftn(F[k][j] * w)) for j in range(3)] for k in (1, 5)}; Dm = Mr[5][0] - Mr[1][0]; Nd = ((Mr[5][1] - Mr[5][2]) - (Mr[1][1] - Mr[1][2])) / 2
                foc = mk & (np.abs(Dm) > 3 * Nd[mk].std())
                if foc.sum() < 0.005 * mk.sum(): foc = mk & (np.abs(Dm) >= np.percentile(np.abs(Dm[mk]), 95))
                use = foc.astype(float)
            sm = np.real(np.fft.ifftn(np.fft.fftn(use) * np.exp(-0.5 * (2. * 2 * np.pi * np.sqrt(np.fft.fftfreq(n)[:, None, None] ** 2 + np.fft.fftfreq(n)[None, :, None] ** 2 + np.fft.fftfreq(n)[None, None, :] ** 2)) ** 2)))
            Q = {k: [np.fft.fftn(np.real(np.fft.ifftn(x)) * sm) for x in F[k]] for k in F}; smax = int(n * smpd / res); sel = (sh >= 1) & (sh <= smax); shs = sh[sel]
            Qs = {k: [x.ravel()[sel] for x in Q[k]] for k in Q}; neff = use.sum() * smpd ** 3 / (res / 2) ** 3; out = []
            for k in (2, 3, 4):
                p1, p5, m = Qs[1][0], Qs[5][0], Qs[k][0]; ab = np.zeros((smax + 1, 2))
                for sidx in range(1, smax + 1):
                    q = shs == sidx; A11 = (p1[q].conj() * p1[q]).real.sum(); A12 = (p1[q].conj() * p5[q]).real.sum(); A22 = (p5[q].conj() * p5[q]).real.sum()
                    y1 = (p1[q].conj() * m[q]).real.sum(); y2 = (p5[q].conj() * m[q]).real.sum(); det = A11 * A22 - A12 * A12
                    ab[sidx] = [(A22 * y1 - A12 * y2) / det, (A11 * y2 - A12 * y1) / det] if det > 0 else [0., 0.]
                a_, b_ = ab[shs, 0], ab[shs, 1]; r = m - a_ * p1 - b_ * p5; nz = lambda kk: np.abs(Qs[kk][1] - Qs[kk][2]) ** 2 / 4
                exp = (nz(k) + a_ ** 2 * nz(1) + b_ ** 2 * nz(5)).sum(); re_ = Qs[k][1] - (a_ * Qs[1][2] + b_ * Qs[5][2]); ro = Qs[k][2] - (a_ * Qs[1][1] + b_ * Qs[5][1])
                rho = float((re_.conj() * ro).real.sum() / np.sqrt((np.abs(re_) ** 2).sum() * (np.abs(ro) ** 2).sum())); pw = np.array([(np.abs(m[shs == q_]) ** 2).sum() for q_ in range(1, smax + 1)])
                bt = ab[1:, 1] / np.where(np.abs(ab[1:].sum(1)) > 1e-9, ab[1:].sum(1), np.nan); okb = np.isfinite(bt); beta = float((bt[okb] * pw[okb]).sum() / pw[okb].sum())
                out.append('bin%d %.2f / %.3f (z %.1f) / b %.2f' % (k, (np.abs(r) ** 2).sum() / exp, rho, rho * np.sqrt(neff), beta))
            print('  %2.0f A %s%s: per-shell mixture resid/noise / rho_res (z) / beta: %s' % (res, var, ' (focus %.3f of mask)' % (use.sum() / mk.sum()) if var == 'SF' else '', ' | '.join(out)))
    sys.stdout.flush()

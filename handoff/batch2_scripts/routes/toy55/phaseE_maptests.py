#!/usr/bin/env python3
"""usage: phaseE_maptests.py NAME:SERIES_DIR:MASK:MSKDIAM[:TRUTHKEY] ..   (phaseD_maptests.py plus the at-node gate and magnitudes)
phaseD's tests: variant A (raw maps in the root mask: step profile, mixture), the per-shell mixture S (soft root mask) and SF (soft
focus mask from the raw P5 - P1, 3 sd of its half-map noise), and where the change sits. Added, per resolution, for S and SF:
- the ends SNR: |P5 - P1|^2 over its half-map noise (|E5 - O5|^2 + |E1 - O1|^2) / 4, summed over shells 1..res; 1 = the ends agree;
- the reproducible ends difference: rho_ends = Re<E5 - E1, O5 - O1> over the norms, and z_ends = rho_ends sqrt(N_eff);
- the non-blend share of each middle bin: the reproducible residual power Re<re, ro> of the cross-fitted per-shell mixture over the
  reproducible ends power Re<E5 - E1, O5 - O1>; ~0 for a blend of the ends whatever the particle count;
- variant A's ends noise ratio, (1 - corr(P1, P5)) over the noise expectation of two maps of one state.
The beta implied by the truth is (t_k - t_1) / (t_5 - t_1), from bins.txt's per-bin open / A / closed fraction or angle mean.
Each series ends with one SUMMARY line, NAME first."""
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
def cdot(x, y): return float((x.conj() * y).real.sum())
def truth(D, key):
    t = []
    for l in open(D + '/bins.txt'):
        if not l.startswith('bin '): continue
        m = re.search(r'%s mean ([+-]?[\d.]+)' % 'a\\d_\\w\\w', l) if key == 'mean' else re.search(r'\b%s ([\d.]+)' % key, l); t.append(float(m.group(1)))
    return np.array(t)
def rows(D): return [int(re.search(r'rows +(\d+)', l).group(1)) for l in open(D + '/bins.txt') if l.startswith('bin ')]
for spec in sys.argv[1:]:
    parts = spec.split(':'); NAME, D, MK, MD = parts[:4]; KEY = parts[4] if len(parts) > 4 and parts[4] else None; F = {}
    for k in range(1, 6):
        F[k] = [np.fft.fftn(read_mrc(f)[0]) for f in ('%s/bin%d.mrc' % (D, k), '%s/halves/bin%d_even.mrc' % (D, k), '%s/halves/bin%d_odd.mrc' % (D, k))]
    smpd = read_mrc('%s/bin1.mrc' % D)[1]; n = F[1][0].shape[0]; mk = fresample(read_mrc(MK)[0], n) > 0.5; sh = np.rint(kr(n)).astype(int).ravel()
    t = truth(D, KEY) if KEY else np.full(5, np.nan); bimp = (t[1:4] - t[0]) / (t[4] - t[0]) if KEY else np.full(3, np.nan); nr = rows(D)
    print('######## %s: %s (box %d, %.3f A/px, root mask %d voxels; rows per bin %s; truth %s per bin %s; truth-implied beta bins 2-4 %s)' % (NAME, D, n, smpd, mk.sum(), ' '.join(map(str, nr)), KEY,
          ' '.join('%.3f' % x for x in t), ' '.join('%.2f' % x for x in bimp)))
    w15 = lp_w(n, smpd, 15.); Dr = np.real(np.fft.ifftn((F[5][0] - F[1][0]) * w15)); Ndr = np.real(np.fft.ifftn(((F[5][1] - F[5][2]) - (F[1][1] - F[1][2])) * w15)) / 2
    focf = (mk & (np.abs(Dr) > 3 * Ndr[mk].std())).sum() / mk.sum(); g = np.arange(n) - n // 2; ZZ, YY, XX = np.meshgrid(g, g, g, indexing='ij')
    def cen(sel): z_, y_, x_ = ZZ[sel].mean() * smpd, YY[sel].mean() * smpd, XX[sel].mean() * smpd; return x_, y_, z_, np.hypot(x_, y_)
    ad = np.abs(Dr); top = mk & (ad >= np.percentile(ad[mk], 99)); pos = mk & (Dr >= np.percentile(Dr[mk], 99.5)); neg = mk & (Dr <= np.percentile(Dr[mk], 0.5))
    ct, cp, cn = cen(top), cen(pos), cen(neg); rr = np.sqrt((XX[top] * smpd) ** 2 + (YY[top] * smpd) ** 2)
    pmr = (Dr[mk] ** 2 * (Dr[mk] > 0)).sum() / max((Dr[mk] ** 2 * (Dr[mk] < 0)).sum(), 1e-30)
    print('  where (raw P5 - P1 at 15 A): focus fraction %.3f of the mask; top 1 %% abs centroid x %.0f y %.0f z %.0f A (radius from z %.0f, median voxel radius %.0f A of mask radius %.0f); + lobe (%.0f, %.0f, %.0f) - lobe (%.0f, %.0f, %.0f), separation %.0f A; + / - power ratio %.2f'
          % (focf, ct[0], ct[1], ct[2], ct[3], np.median(rr), float(MD) / 2, cp[0], cp[1], cp[2], cn[0], cn[1], cn[2], np.sqrt(sum((cp[i] - cn[i]) ** 2 for i in range(3))), pmr))
    S = {}
    for res in (15., 10.):
        w = lp_w(n, smpd, res)
        M = {k: [np.real(np.fft.ifftn(x * w)) for x in F[k]] for k in F}; V = {k: [x[mk] for x in M[k]] for k in M}; neff = mk.sum() * smpd ** 3 / (res / 2) ** 3
        eo = {k: 1 - corr(V[k][1], V[k][2]) for k in V}; d = [1 - corr(V[k][0], V[k + 1][0]) for k in range(1, 5)]; rn = [d[k - 1] / ((eo[k] + eo[k + 1]) / 4) for k in range(1, 5)]
        endsA = (1 - corr(V[1][0], V[5][0])) / ((eo[1] + eo[5]) / 4); P1, P5 = V[1][0], V[5][0]; nv = lambda k: np.var(V[k][1] - V[k][2]) / 4; out = []
        for k in (2, 3, 4):
            A = np.column_stack([P1, P5, np.ones(len(P1))]); a, b, c = np.linalg.lstsq(A, V[k][0], rcond=None)[0]; rr = V[k][0] - A @ np.array([a, b, c])
            exp = nv(k) + a * a * nv(1) + b * b * nv(5); re_ = V[k][1] - (a * V[1][2] + b * V[5][2] + c); ro = V[k][2] - (a * V[1][1] + b * V[5][1] + c); rho = corr(re_, ro)
            beta = np.linalg.lstsq(np.column_stack([P5 - P1, np.ones(len(P1))]), V[k][0] - P1, rcond=None)[0][0]
            out.append('bin%d %.2f / %.3f (z %.1f) / b %.2f' % (k, np.var(rr) / exp, rho, rho * np.sqrt(neff), beta))
        print('  %2.0f A A: steps %s, dominant share %.2f, ends noise ratio %.1f | mixture resid/noise / rho_res (z) / beta: %s' % (res, ' '.join('%.1f' % x for x in rn), max(rn) / sum(rn), endsA, ' | '.join(out)))
        S[res, 'step'] = max(rn) / sum(rn); S[res, 'endsA'] = endsA; S[res, 'steps'] = rn
        for var in ('S', 'SF'):
            if var == 'S': use = mk.astype(float)
            else:
                Mr = {k: [np.real(np.fft.ifftn(F[k][j] * w)) for j in range(3)] for k in (1, 5)}; Dm = Mr[5][0] - Mr[1][0]; Nd = ((Mr[5][1] - Mr[5][2]) - (Mr[1][1] - Mr[1][2])) / 2
                foc = mk & (np.abs(Dm) > 3 * Nd[mk].std())
                if foc.sum() < 0.005 * mk.sum(): foc = mk & (np.abs(Dm) >= np.percentile(np.abs(Dm[mk]), 95))
                use = foc.astype(float)
            sm = np.real(np.fft.ifftn(np.fft.fftn(use) * np.exp(-0.5 * (2. * 2 * np.pi * np.sqrt(np.fft.fftfreq(n)[:, None, None] ** 2 + np.fft.fftfreq(n)[None, :, None] ** 2 + np.fft.fftfreq(n)[None, None, :] ** 2)) ** 2)))
            Q = {k: [np.fft.fftn(np.real(np.fft.ifftn(x)) * sm) for x in F[k]] for k in F}; smax = int(n * smpd / res); sel = (sh >= 1) & (sh <= smax); shs = sh[sel]
            Qs = {k: [x.ravel()[sel] for x in Q[k]] for k in Q}; neff = use.sum() * smpd ** 3 / (res / 2) ** 3; out = []; R = []; NB = []
            de, do = Qs[5][1] - Qs[1][1], Qs[5][2] - Qs[1][2]; nd = ((np.abs(Qs[5][1] - Qs[5][2]) ** 2 + np.abs(Qs[1][1] - Qs[1][2]) ** 2) / 4).sum()
            esnr = (np.abs(Qs[5][0] - Qs[1][0]) ** 2).sum() / nd; rep = cdot(de, do); rhoe = rep / np.sqrt((np.abs(de) ** 2).sum() * (np.abs(do) ** 2).sum())
            for k in (2, 3, 4):
                p1, p5, m = Qs[1][0], Qs[5][0], Qs[k][0]; ab = np.zeros((smax + 1, 2))
                for sidx in range(1, smax + 1):
                    q = shs == sidx; A11 = (p1[q].conj() * p1[q]).real.sum(); A12 = (p1[q].conj() * p5[q]).real.sum(); A22 = (p5[q].conj() * p5[q]).real.sum()
                    y1 = (p1[q].conj() * m[q]).real.sum(); y2 = (p5[q].conj() * m[q]).real.sum(); det = A11 * A22 - A12 * A12
                    ab[sidx] = [(A22 * y1 - A12 * y2) / det, (A11 * y2 - A12 * y1) / det] if det > 0 else [0., 0.]
                a_, b_ = ab[shs, 0], ab[shs, 1]; r = m - a_ * p1 - b_ * p5; nz = lambda kk: np.abs(Qs[kk][1] - Qs[kk][2]) ** 2 / 4
                exp = (nz(k) + a_ ** 2 * nz(1) + b_ ** 2 * nz(5)).sum(); re_ = Qs[k][1] - (a_ * Qs[1][2] + b_ * Qs[5][2]); ro = Qs[k][2] - (a_ * Qs[1][1] + b_ * Qs[5][1])
                rho = cdot(re_, ro) / np.sqrt((np.abs(re_) ** 2).sum() * (np.abs(ro) ** 2).sum()); pw = np.array([(np.abs(m[shs == q_]) ** 2).sum() for q_ in range(1, smax + 1)])
                bt = ab[1:, 1] / np.where(np.abs(ab[1:].sum(1)) > 1e-9, ab[1:].sum(1), np.nan); okb = np.isfinite(bt); beta = float((bt[okb] * pw[okb]).sum() / pw[okb].sum()) if okb.any() else np.nan
                nb = cdot(re_, ro) / rep if rep > 0 else np.nan; R.append((np.abs(r) ** 2).sum() / exp); NB.append(nb)
                out.append('bin%d %.2f / %.3f (z %.1f) / b %.2f / nb %.3f' % (k, R[-1], rho, rho * np.sqrt(neff), beta, nb))
            print('  %2.0f A %s%s: ends SNR %.2f, rho_ends %.3f (z %.1f) | per-shell mixture resid/noise / rho_res (z) / beta / non-blend share: %s' % (res, var,
                  ' (focus %.3f of mask)' % (use.sum() / mk.sum()) if var == 'SF' else '', esnr, rhoe, rhoe * np.sqrt(neff), ' | '.join(out)))
            S[res, var] = (esnr, rhoe * np.sqrt(neff), R, NB)
    def fmt(res):
        e, z, R, NB = S[res, 'S']; e2, z2, R2, NB2 = S[res, 'SF']
        return '%2.0f A: ends SNR %.2f z %.1f (A ratio %.1f) steps %s share %.2f | S %s nb %s | SF %s nb %s' % (res, e, z, S[res, 'endsA'], '/'.join('%.1f' % x for x in S[res, 'steps']), S[res, 'step'],
               '/'.join('%.2f' % x for x in R), '/'.join('%.3f' % x for x in NB), '/'.join('%.2f' % x for x in R2), '/'.join('%.3f' % x for x in NB2))
    print('SUMMARY %s | min rows/bin %d | focus %.2f, +/- %.2f | %s | %s' % (NAME, min(nr), focf, pmr, fmt(15.), fmt(10.)))
    sys.stdout.flush()

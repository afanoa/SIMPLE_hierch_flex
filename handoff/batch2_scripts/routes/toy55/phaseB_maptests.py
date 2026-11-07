#!/usr/bin/env python3
"""usage: phaseB_maptests.py SERIES_DIR:MASK:MSKDIAM ..
Map-level tests of a 5-bin series (bin1..5.mrc, halves/binK_even|odd.mrc) inside the tree's root MASK (Fourier-resampled to the map
box), low-passed at 15 and 10 A.
(1) Step profile: d_k = 1 - corr(M_k, M_k+1); its noise ratio to [(1 - corr(E_k, O_k)) + (1 - corr(E_k+1, O_k+1))] / 4, the distance two
    full maps of one state would show; the dominant step's share of the summed ratios.
(2) Mixture test for k = 2, 3, 4: M_k = a P1 + b P5 + c by least squares (P1 / P5 = bins 1 / 5). Reported:
  - a, b, a + b;
  - residual variance / the noise expectation, var(E_k - O_k) / 4 + a^2 var(E1 - O1) / 4 + b^2 var(E5 - O5) / 4;
  - the reproducible residual rho_res = corr(E_k - (a P1_o + b P5_o + c), O_k - (a P1_e + b P5_e + c)), and z = rho_res sqrt(N_eff)
    with N_eff = mask volume / (res / 2)^3;
  - the reproducible residual power in the moving region (top 10 % of abs(P5 - P1)), as an enrichment over its 10 % volume;
  - beta of the constrained fit M_k = (1 - beta) P1 + beta P5 + c."""
import sys, numpy as np
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
def lowpass(v, smpd, res):
    n = v.shape[0]; f = np.fft.fftfreq(n); r = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
    fc = smpd / res; w = np.clip((1.1 * fc - r) / (0.2 * fc), 0., 1.); w = 0.5 - 0.5 * np.cos(np.pi * w); return np.real(np.fft.ifftn(np.fft.fftn(v) * w))
def corr(a, b): a = a - a.mean(); b = b - b.mean(); return float(a @ b / np.sqrt((a @ a) * (b @ b)))
for spec in sys.argv[1:]:
    D, MK, MD = spec.split(':'); maps = {}
    for k in range(1, 6):
        maps[k] = [read_mrc('%s/bin%d.mrc' % (D, k))[0], read_mrc('%s/halves/bin%d_even.mrc' % (D, k))[0], read_mrc('%s/halves/bin%d_odd.mrc' % (D, k))[0]]
    smpd = read_mrc('%s/bin1.mrc' % D)[1]; box = maps[1][0].shape[0]; mk = fresample(read_mrc(MK)[0], box) > 0.5
    print('######## %s (box %d, %.3f A/px, mask %d voxels)' % (D, box, smpd, mk.sum()))
    for res in (15., 10.):
        V = {k: [lowpass(x, smpd, res)[mk] for x in maps[k]] for k in maps}; neff = mk.sum() * smpd ** 3 / (res / 2) ** 3
        eo = {k: 1 - corr(V[k][1], V[k][2]) for k in V}; d = [1 - corr(V[k][0], V[k + 1][0]) for k in range(1, 5)]
        rn = [d[k - 1] / ((eo[k] + eo[k + 1]) / 4) for k in range(1, 5)]
        print('  %2.0f A step profile: 1-corr %s | noise ratio %s | dominant share %.2f | ends 1-corr %.4f (noise ratio %.1f)' % (res, ' '.join('%.4f' % x for x in d),
              ' '.join('%.1f' % x for x in rn), max(rn) / sum(rn), 1 - corr(V[1][0], V[5][0]), (1 - corr(V[1][0], V[5][0])) / ((eo[1] + eo[5]) / 4)))
        P1, P5 = V[1][0], V[5][0]; mov = np.abs(P5 - P1) >= np.percentile(np.abs(P5 - P1), 90); nv = lambda k: np.var(V[k][1] - V[k][2]) / 4
        for k in (2, 3, 4):
            A = np.column_stack([P1, P5, np.ones(len(P1))]); a, b, c = np.linalg.lstsq(A, V[k][0], rcond=None)[0]; rr = V[k][0] - A @ np.array([a, b, c])
            exp = nv(k) + a * a * nv(1) + b * b * nv(5); re = V[k][1] - (a * V[1][2] + b * V[5][2] + c); ro = V[k][2] - (a * V[1][1] + b * V[5][1] + c)
            rho = corr(re, ro); pw = (re - re.mean()) * (ro - ro.mean()); enr = (pw[mov].sum() / pw.sum()) / mov.mean() if pw.sum() > 0 else float('nan')
            beta = np.linalg.lstsq(np.column_stack([P5 - P1, np.ones(len(P1))]), V[k][0] - P1, rcond=None)[0][0]
            print('  %2.0f A mixture bin%d: a %.3f b %.3f a+b %.3f | resid var / noise %.2f | reproducible resid rho %.3f, z %.1f | moving-region enrichment %.2f | constrained beta %.3f'
                  % (res, k, a, b, a + b, np.var(rr) / exp, rho, rho * np.sqrt(neff), enr, beta))
    sys.stdout.flush()

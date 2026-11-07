#!/usr/bin/env python3
"""usage: phaseF_rel.py NAME:SERIES_DIR:MASK ..   (how big and how local a 5-bin series' end-to-end change is)
In phaseE_maptests.py's S frame (soft root mask, shells 1..15 A and 1..10 A), with D = P5 - P1 and the mean map Pm over the 5 bins:
- rel: the reproducible end-to-end power Re<E5 - E1, O5 - O1> over the mean map's, Re<Em, Om>;
- cos: the signal cosine of D with Pm, Re<E5 - E1, Om> over sqrt(both reproducible powers); 1 = D is a rescaled copy of the map
  (brightness, junk dilution), 0 = D is unrelated to it;
- the local change: D less its per-shell projection on Pm, c(s) = Re<D, Pm> / |Pm|^2 (fitted on the full maps, applied to each half):
  rel_loc (its reproducible power over the mean map's) and z_loc (the split-half correlation of the local change times sqrt(N_eff));
- hmin: the worst bin's half-map correlation, Re<E_k, O_k> over the norms (how well each bin's own map is determined);
- res: the worst bin's FSC 0.143 resolution in the soft mask (the first shell below 0.143), A.
One RELF line per series."""
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
def kr(n): f = np.fft.fftfreq(n) * n; return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def cdot(x, y): return float((x.conj() * y).real.sum())
for spec in sys.argv[1:]:
    NAME, D, MK = spec.split(':')[:3]; R = {}
    V = {k: [read_mrc(f)[0] for f in ('%s/bin%d.mrc' % (D, k), '%s/halves/bin%d_even.mrc' % (D, k), '%s/halves/bin%d_odd.mrc' % (D, k))] for k in range(1, 6)}
    smpd = read_mrc('%s/bin1.mrc' % D)[1]; n = V[1][0].shape[0]; mk = (fresample(read_mrc(MK)[0], n) > 0.5).astype(float); sh = np.rint(kr(n)).astype(int).ravel()
    sm = np.real(np.fft.ifftn(np.fft.fftn(mk) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2)))
    Q = {k: [np.fft.fftn(x * sm).ravel() for x in V[k]] for k in V}; ns = n // 2; rb = []
    for k in Q:
        ee = np.bincount(sh, np.abs(Q[k][1]) ** 2, ns + 1); oo = np.bincount(sh, np.abs(Q[k][2]) ** 2, ns + 1); eo = np.bincount(sh, (Q[k][1] * Q[k][2].conj()).real, ns + 1)
        fsc = eo[1:] / np.sqrt(np.maximum(ee[1:] * oo[1:], 1e-30)); lo = np.nonzero(fsc < 0.143)[0]; rb.append(n * smpd / (lo[0] + 1) if len(lo) else 2 * smpd)
    for res in (15., 10.):
        smax = int(n * smpd / res); sel = (sh >= 1) & (sh <= smax); shs = sh[sel]; neff = mk.sum() * smpd ** 3 / (res / 2) ** 3
        q = {k: [x[sel] for x in Q[k]] for k in Q}; pm = [np.mean([q[k][j] for k in q], 0) for j in range(3)]
        d = [q[5][j] - q[1][j] for j in range(3)]; sig_m = cdot(pm[1], pm[2]); sig_e = cdot(d[1], d[2])
        cs = cdot(d[1], pm[2]) / np.sqrt(max(sig_e, 1e-30) * max(sig_m, 1e-30)) if sig_e > 0 else np.nan
        num = np.bincount(shs, (d[0].conj() * pm[0]).real, smax + 1); den = np.bincount(shs, np.abs(pm[0]) ** 2, smax + 1); c = (num / np.maximum(den, 1e-30))[shs]
        le, lo = d[1] - c * pm[1], d[2] - c * pm[2]; sig_l = cdot(le, lo); rl = sig_l / np.sqrt(cdot(le, le) * cdot(lo, lo))
        hmin = min(cdot(q[k][1], q[k][2]) / np.sqrt(cdot(q[k][1], q[k][1]) * cdot(q[k][2], q[k][2])) for k in q)
        R[res] = (sig_e / sig_m, cs, sig_l / sig_m, rl * np.sqrt(neff), hmin)
    print('RELF %s | worst-bin FSC 0.143 %.1f A | ' % (NAME, max(rb)) + ' | '.join('%2.0f A: rel %.4f cos %+.3f rel_loc %.4f z_loc %.1f hmin %.3f' % ((res,) + R[res]) for res in (15., 10.)))
    sys.stdout.flush()

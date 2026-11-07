#!/usr/bin/env python3
"""usage: leaf_stats.py TREE_DIR DUMP_NPZ MAPLIST OUT
Per leaf of a hetero_tree run (MAPLIST: fa_maplist_final.sh lines 'name kind family particles map even odd', leaves and junk units):
- the particles labelled to it (hetero_tree_labels.txt index -> hetero_tree_leaves.txt line);
- quality: its corr against all other root rows at matched defocus (dfx deciles, root stage-B corr from DUMP_NPZ; 0.5 = equal),
  its median dfx and top-view share;
- FSC 0.143 of its halves in the soft root mask (TREE_DIR/node_r/mask.mrc), A.
One LEAF line per map, written to OUT."""
import sys, numpy as np
T, DP, ML, OUT = sys.argv[1:5]; T = T.rstrip('/')
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
D = np.load(DP); corr, dfx, e2 = D['corr'], D['dfx'], D['e2']; ok = np.isfinite(dfx)
dd = np.full(len(dfx), -1); dd[ok] = np.searchsorted(np.quantile(dfx[ok], np.linspace(0.1, 0.9, 9)), dfx[ok])
L = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int); LV = [l.split()[0] for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]; idx = {nm: k + 1 for k, nm in enumerate(LV)}
assert len(L) == len(corr), 'labels %d rows, dump %d' % (len(L), len(corr))
def auc(x, y):
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 20 or len(y) < 20: return np.nan
    r = np.argsort(np.argsort(np.r_[x, y])); return (r[:len(x)].sum() - len(x) * (len(x) - 1) / 2.) / (len(x) * len(y))
mk0 = read_mrc(T + '/node_r/mask.mrc')[0]; inroot = L > 0; cache = {}
with open(OUT, 'w') as fo:
    for w in [l.split() for l in open(ML) if l.strip() and not l.startswith('#')]:
        nm, kind = w[0], w[1]; k = idx.get(nm, -1); m = L == k
        if kind not in ('leaf', 'junk'): continue
        q = np.nanmean([auc(corr[m & (dd == j)], corr[inroot & ~m & (dd == j)]) for j in range(10)]) if m.sum() >= 20 else np.nan
        res = np.nan
        if len(w) >= 7:
            (e, smpd), (o, _) = read_mrc(w[5]), read_mrc(w[6]); n = e.shape[0]
            if n not in cache:
                mk = (fresample(mk0, n) > 0.5).astype(float); cache[n] = (np.real(np.fft.ifftn(np.fft.fftn(mk) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2))), np.rint(kr(n)).astype(int).ravel())
            sm, sh = cache[n]; E, O = np.fft.fftn(e * sm).ravel(), np.fft.fftn(o * sm).ravel(); ns = n // 2
            eo = np.bincount(sh, (E * O.conj()).real, ns + 1); ee = np.bincount(sh, np.abs(E) ** 2, ns + 1); oo = np.bincount(sh, np.abs(O) ** 2, ns + 1)
            fsc = eo[1:] / np.sqrt(np.maximum(ee[1:] * oo[1:], 1e-30)); lo = np.nonzero(fsc < 0.143)[0]; res = n * smpd / (lo[0] + 1) if len(lo) else 2 * smpd
        fo.write('LEAF %s %s | rows %d | qAUC %.3f | median dfx %.3f | top-view share %.3f | FSC0.143 %.2f A\n' % (nm, kind, m.sum(), q, np.nanmedian(dfx[m]) if m.any() else np.nan,
                 np.mean(np.abs(np.cos(np.radians(e2[m]))) > 0.8) if m.any() else np.nan, res)); fo.flush()

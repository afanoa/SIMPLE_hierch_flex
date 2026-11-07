#!/usr/bin/env python3
"""usage: paths.py TRUTH OUT NAME:TREE_DIR:TABLE:FLEX_EMB:FFLEX_EMB ..
Do a tree's states, or a continuous score, give the motion extremes AND the way between them? (IgM arm angles; genuine particles;
each angle standardised over all genuine particles.)
(1) State paths per motion. The ends are the states (>= 2,000 genuine) with the lowest and the highest mean angle; states with a
    junk share > 0.5 or < 200 genuine members are left out.
  (A) PAGA-like latent graph: a 15-NN graph of genuine particles in the root flex latent (z standardised). Pair connectivity = the
      edges between two states / the 30 n_a n_b / N expected under random mixing. Pairs above 1 are edges, and the path is the
      shortest one on 1 / connectivity.
  (B) Map graph: stage-B leaf maps (a unit = the population-weighted mean of its leaves' maps), low-passed at 15 A inside the root
      mask. Weight 1 - correlation; the path runs along the minimum spanning tree.
  Per path: the number of states, the Spearman of position vs the states' mean angle, and the largest step between neighbours (sd).
(2) Continuous scores: per kept site, bins of the pole-1 posterior with edges at genuine quantiles (5 equal bins; and 10k end bins +
    3 equal middle bins). Ties are broken at random, and junk is counted in the same score intervals.
  - Per bin: the mean / sd of the site's best-matching angle (oriented to rise with the score) and the junk share.
  - Per site: the Spearman of bin index vs mean, the end-to-end span (sd), its ratio to the oracle's span at the end bins' genuine
    sizes, and the largest / smallest step.
  - The same for the root flex / fflex component that best matches each angle.
(3) Coverage: per angle, the best abs(Spearman) of a kept site and of any candidate posterior (node site_post_*.txt and fsite_*)."""
import sys, os, glob, numpy as np
from scipy.spatial import cKDTree
TF, OUT = sys.argv[1:3]; SPECS = [a.split(':') for a in sys.argv[3:]]; RNG = np.random.default_rng(7)
hd = [l.split()[1:] for l in open(TF) if l.startswith('#')][0]; T = np.genfromtxt(TF, comments='#'); tr = T[:, 0].astype(int) - 1; NM = tr.max() + 1
ang = hd[3:]; NA = len(ang); X = np.full((NM, NA), np.nan); X[tr] = T[:, 3:3 + NA]
junk = np.zeros(NM, bool); junk[tr] = T[:, hd.index('junk')] > 0; G = ~junk & ~np.isnan(X[:, 0])
Z = np.nan_to_num((X - X[G].mean(0)) / X[G].std(0))
EH = [np.cumsum(np.sort(Z[G, c])[::-1]) / np.arange(1, G.sum() + 1) for c in range(NA)]       # oracle mean of the n highest
EL = [np.cumsum(np.sort(-Z[G, c])[::-1]) / np.arange(1, G.sum() + 1) for c in range(NA)]      # ... of the n lowest (as a positive sd)
def rank(v): r = np.empty(len(v)); r[np.argsort(v, kind='stable')] = np.arange(len(v)); return r
def spear(x, y): return float(np.corrcoef(rank(np.asarray(x, float)), rank(np.asarray(y, float)))[0, 1]) if len(x) > 2 else float(np.sign(y[-1] - y[0]))
def read_mrc(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = np.frombuffer(raw[:1024], np.float32)
    nx, ny, nz, mode, nsym = int(h[0]), int(h[1]), int(h[2]), int(h[3]), int(h[23])
    v = np.frombuffer(raw[1024 + nsym:], {0: np.int8, 1: np.int16, 2: np.float32, 6: np.uint16}[mode], count=nx * ny * nz).astype(np.float64)
    return v.reshape(nz, ny, nx), float(hf[10]) / nx
def fresample(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2
    if n < m: s = c - n // 2; Gf = F[s:s + n, s:s + n, s:s + n]
    else: Gf = np.zeros((n, n, n), complex); s = n // 2 - c; Gf[s:s + m, s:s + m, s:s + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(Gf))) * (n / m) ** 3
def lowpass(v, smpd, res):
    n = v.shape[0]; f = np.fft.fftfreq(n); r = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
    fc = smpd / res; w = np.clip((1.1 * fc - r) / (0.2 * fc), 0., 1.); w = 0.5 - 0.5 * np.cos(np.pi * w)
    return np.real(np.fft.ifftn(np.fft.fftn(v) * w))
def read_emb(f):
    b = open(f, 'rb').read(); ver, n, nc, box = (int(v) for v in np.frombuffer(b[8:24], '<i4'))
    assert b[:8] == b'SIMPLFXC' and ver == 5, 'flex_pca embedding version 5 expected'
    o = 28; pinds = np.frombuffer(b[o:o + 4 * n], '<i4').astype(int); o += 4 * n
    return pinds - 1, np.frombuffer(b[o:o + 8 * n * nc], '<f8').reshape(nc, n).T.copy()
def readcol(f): return np.array(open(f).read().split(), float)
def dijkstra(W, a, b):
    m = len(W); d = np.full(m, np.inf); d[a] = 0.; prev = -np.ones(m, int); done = np.zeros(m, bool)
    while True:
        dm = np.where(done, np.inf, d); u = int(np.argmin(dm))
        if not np.isfinite(dm[u]) or u == b: break
        done[u] = True; nd = d[u] + W[u]; up = (~done) & (nd < d); d[up] = nd[up]; prev[up] = u
    if not np.isfinite(d[b]): return None
    p = [b]
    while p[-1] != a: p.append(prev[p[-1]])
    return p[::-1]
def mst_path(Wt, a, b):
    m = len(Wt); intree = np.zeros(m, bool); intree[0] = True; best = Wt[0].copy(); par = np.zeros(m, int); adj = {i: [] for i in range(m)}
    for _ in range(m - 1):
        u = int(np.argmin(np.where(intree, np.inf, best))); intree[u] = True; adj[u].append(par[u]); adj[par[u]].append(u)
        up = (~intree) & (Wt[u] < best); best[up] = Wt[u][up]; par[up] = u
    prev = {a: None}; q = [a]
    while q:
        u = q.pop(0)
        for v in adj[u]:
            if v not in prev: prev[v] = u; q.append(v)
    p = [b]
    while p[-1] != a: p.append(prev[p[-1]])
    return p[::-1]
def bins_stats(score, rows, c, edges_kind):
    gl = rows & G; s = score + RNG.uniform(0., 1e-9, NM); sg = np.sort(s[gl]); ng = len(sg)
    if edges_kind == 'q5': cuts = [sg[int(ng * q) - 1] for q in (0.2, 0.4, 0.6, 0.8)]
    else:
        if ng < 25000: return None
        mid = (ng - 20000) // 3; cuts = [sg[10000 - 1], sg[10000 + mid - 1], sg[10000 + 2 * mid - 1], sg[ng - 10000 - 1]]
    b = np.searchsorted(np.array(cuts), s, side='left'); out = []
    for k in range(5):
        m = rows & (b == k); g = m & G
        out.append((g.sum(), Z[g, c].mean() if g.any() else np.nan, Z[g, c].std() if g.any() else np.nan, (m & junk).sum() / max(m.sum(), 1)))
    return out
def summarise(bs, c, sgn):
    mu = np.array([x[1] * sgn for x in bs]); n0, n4 = bs[0][0], bs[-1][0]; span = mu[-1] - mu[0]
    ospan = (EH[c] if sgn > 0 else EL[c])[n4 - 1] + (EL[c] if sgn > 0 else EH[c])[n0 - 1]; st = np.diff(mu)
    return spear(np.arange(5), mu), span, span / ospan, st.max(), st.min(), mu
P1, P2, P3 = open(OUT + '_paths.tsv', 'w'), open(OUT + '_bins.tsv', 'w'), open(OUT + '_cover.tsv', 'w')
P1.write('tree\tstates\tgraph\tangle\tn_states\tspearman\tmax_step\tend_low\tend_high\tpath\n')
P2.write('tree\tsource\tangle\tsign\tabs_rho\tbins\tspearman\tspan\tspan_ratio\tmax_step\tmin_step\tbin_means\tbin_sd\tbin_genuine\tbin_junk\n')
P3.write('tree\tangle\tbest_kept\tbest_kept_site\tbest_candidate\tbest_candidate_file\tcandidates_ge_0.5\tcandidates\n')
SUM = []
for sp in SPECS:
    nm, D, tab, fe, ffe = sp; L = np.loadtxt(D + '/hetero_tree_labels.txt', dtype=int)[:NM]
    LV = [l.split() for l in open(D + '/hetero_tree_leaves.txt') if l.strip()]
    U = np.loadtxt(D + '/hetero_tree_merged.txt', dtype=int)[:NM] if os.path.exists(D + '/hetero_tree_merged.txt') else L.copy()
    mask, msp = read_mrc(D + '/node_r/mask.mrc'); leafmap = {}
    def lmap(k):
        if k not in leafmap: v, smpd = read_mrc(LV[k - 1][4]); leafmap[k] = (v, smpd)
        return leafmap[k]
    lu = {k: int(np.bincount(U[L == k]).argmax()) for k in np.unique(L[L > 0])}
    R, z = read_emb(fe); g = G[R]; Rg = R[g]; Zl = z[g]; Zl = (Zl - Zl.mean(0)) / Zl.std(0)
    nb = cKDTree(Zl).query(Zl, k=16, workers=2)[1][:, 1:]; print('%s: 15-NN graph of %d genuine particles' % (nm, len(Rg)), flush=True)
    for kind, S in (('leaves', L), ('units', U)):
        ids = [s for s in np.unique(S[S > 0]) if (S == s).sum() and junk[S == s].mean() <= 0.5 and ((S == s) & G).sum() >= 200]
        ix = -np.ones(S.max() + 1, int); ix[ids] = np.arange(len(ids)); m = len(ids)
        ng = np.array([((S == s) & G).sum() for s in ids]); mu = np.array([Z[(S == s) & G].mean(0) for s in ids])
        names = [LV[s - 1][0] if kind == 'leaves' else 'u%d' % s for s in ids]
        la = ix[S[Rg]]; lb = la[nb]; ok = (la[:, None] >= 0) & (lb >= 0); C = np.zeros((m, m)); np.add.at(C, (np.broadcast_to(la[:, None], lb.shape)[ok], lb[ok]), 1)
        Cs = C + C.T; npg = np.bincount(la[la >= 0], minlength=m).astype(float); N = npg.sum(); conn = Cs / np.maximum(30. * np.outer(npg, npg) / N, 1e-12)
        WA = np.where((conn > 1.) & ~np.eye(m, dtype=bool), 1. / np.maximum(conn, 1e-12), np.inf)
        WA0 = np.where((conn > 0.) & ~np.eye(m, dtype=bool), 1. / np.maximum(conn, 1e-12), np.inf)   # every observed link (A0)
        vecs = []
        for s in ids:
            mem = [k for k in lu if (lu[k] == s if kind == 'units' else k == s)]; wts = np.array([(L == k).sum() for k in mem], float)
            v0, smpd = lmap(mem[0]); acc = np.zeros_like(v0)
            for k, w in zip(mem, wts): acc += w * lmap(k)[0]
            vecs.append((acc / wts.sum(), smpd))
        box, smpd = vecs[0][0].shape[0], vecs[0][1]; mk = fresample(mask, box) > 0.5
        V = np.array([lowpass(v, smpd, 15.)[mk] for v, _ in vecs]); V = (V - V.mean(1, keepdims=True)); V /= np.linalg.norm(V, axis=1, keepdims=True)
        WB = 1. - V @ V.T; np.fill_diagonal(WB, np.inf)
        big = np.where(ng >= 2000)[0]; res = {'A': [], 'A0': [], 'B': []}
        for c in range(NA):
            lo, hi = big[np.argmin(mu[big, c])], big[np.argmax(mu[big, c])]
            for gname, path in (('A', dijkstra(WA, lo, hi)), ('A0', dijkstra(WA0, lo, hi)), ('B', mst_path(WB, lo, hi))):
                if path is None: P1.write('%s\t%s\t%s\t%s\t0\tnan\tnan\t%s\t%s\tdisconnected\n' % (nm, kind, gname, ang[c], names[lo], names[hi])); res[gname].append((0, np.nan, np.nan)); continue
                y = mu[path, c]; r = spear(np.arange(len(path)), y); stp = np.abs(np.diff(y)).max()
                res[gname].append((len(path), r, stp))
                P1.write('%s\t%s\t%s\t%s\t%d\t%.3f\t%.3f\t%s (%.2f)\t%s (%.2f)\t%s\n' % (nm, kind, gname, ang[c], len(path), r, stp, names[lo], mu[lo, c], names[hi], mu[hi, c], ' > '.join('%s(%.2f)' % (names[q], mu[q, c]) for q in path)))
        for gname in ('A', 'A0', 'B'):
            ns = [x[0] for x in res[gname]]; rs = [x[1] for x in res[gname] if x[0]]; sts = [x[2] for x in res[gname] if x[0]]
            ne = {'A': '%d' % int((WA < np.inf).sum() / 2), 'A0': '%d' % int((WA0 < np.inf).sum() / 2), 'B': 'MST'}[gname]
            SUM.append('%-10s %-6s graph %-2s (%d states, %s edges): median states %.1f, median Spearman %.2f (%d of %d motions >= 0.8), median max step %.2f sd, %d disconnected'
                       % (nm, kind, gname, m, ne, np.median([n for n in ns if n]) if rs else np.nan, np.median(rs) if rs else np.nan, sum(r >= 0.8 for r in rs), NA, np.median(sts) if sts else np.nan, sum(n == 0 for n in ns)))
    print('%s: state paths done' % nm, flush=True)
    kept = []
    tp = tab if tab.startswith('/') else D + '/' + tab
    for w in [l.split() for l in open(tp) if l.strip() and not l.startswith('#')]:
        Q = np.loadtxt(w[3], ndmin=2); Q = Q if Q.shape[0] == NM else Q.T; sc = Q[:, 0]; rows = sc >= 0; gl = rows & G
        rh = np.array([spear(sc[gl], Z[gl, c]) for c in range(NA)]); c = int(np.argmax(np.abs(rh))); sg = 1. if rh[c] > 0 else -1.; kept.append((w[0], np.abs(rh)))
        for bk in ('q5', 'end10k'):
            bs = bins_stats(sc, rows, c, bk)
            if bs is None: continue
            r, span, ratio, mx, mn, mus = summarise(bs, c, sg)
            P2.write('%s\tsite %s (%d genuine)\t%s\t%+d\t%.3f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n' % (nm, w[0], gl.sum(), ang[c], sg, abs(rh[c]), bk, r, span, ratio, mx, mn,
                     ' '.join('%.2f' % x for x in mus), ' '.join('%.2f' % x[2] for x in bs), ' '.join('%d' % x[0] for x in bs), ' '.join('%.3f' % x[3] for x in bs)))
    allr = np.ones(NM, bool)
    for lname, ef in (('flex', fe), ('fflex', ffe)):
        R2, z2 = read_emb(ef); full = np.full((NM, z2.shape[1]), np.nan); full[R2] = z2; ok = ~np.isnan(full[:, 0]) & G
        rz = np.array([[spear(full[ok, j], Z[ok, c]) for j in range(z2.shape[1])] for c in range(NA)])
        for c in range(NA):
            j = int(np.argmax(np.abs(rz[c]))); sg = 1. if rz[c, j] > 0 else -1.; sc = np.nan_to_num(full[:, j], nan=-1e30); rows = ~np.isnan(full[:, j])
            for bk in ('q5', 'end10k'):
                bs = bins_stats(sc, rows, c, bk)
                if bs is None: continue
                r, span, ratio, mx, mn, mus = summarise(bs, c, sg)
                P2.write('%s\t%s z%d\t%s\t%+d\t%.3f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n' % (nm, lname, j + 1, ang[c], sg, abs(rz[c, j]), bk, r, span, ratio, mx, mn,
                         ' '.join('%.2f' % x for x in mus), ' '.join('%.2f' % x[2] for x in bs), ' '.join('%d' % x[0] for x in bs), ' '.join('%.3f' % x[3] for x in bs)))
    print('%s: bins done' % nm, flush=True)
    cand = sorted(glob.glob(D + '/node_*/site/site_post_*.txt') + glob.glob(D + '/node_*/fsite_*/site_post_*.txt')); best = np.zeros(NA); bestf = [''] * NA; ge = np.zeros(NA, int); cache = {}
    for f in cand:
        try: sc = readcol(f)
        except Exception: continue
        if len(sc) != NM: continue
        rows = sc >= 0; gl = rows & G
        if gl.sum() < 1000: continue
        key = (f.split('/node_')[1].split('/')[0], int(gl.sum()))
        if key not in cache: cache = {key: np.array([rank(Z[gl, c]) for c in range(NA)])}
        rs = rank(sc[gl]); ar = np.abs([np.corrcoef(rs, cache[key][c])[0, 1] for c in range(NA)])
        up = ar > best; best[up] = ar[up]; ge += ar >= 0.5
        for c in np.where(up)[0]: bestf[c] = f.split('/1_hetero_tree/')[1]
    for c in range(NA):
        kb = max(kept, key=lambda x: x[1][c]) if kept else ('-', np.zeros(NA))
        P3.write('%s\t%s\t%.3f\t%s\t%.3f\t%s\t%d\t%d\n' % (nm, ang[c], kb[1][c], kb[0], best[c], bestf[c], ge[c], len(cand)))
    print('%s: coverage done over %d candidate posteriors' % (nm, len(cand)), flush=True)
P1.close(); P2.close(); P3.close()
print('\n'.join(SUM))

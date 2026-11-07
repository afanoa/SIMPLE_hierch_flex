#!/usr/bin/env python3
"""Phase-0 scorecard of a state tree (rflex13 / rflex13s run dir) against per-site truth (site_truth.py), per SITE:
  - NMI(site state; leaf) and, for continuous sites, eta2 of the site variable over the leaves (non-junk rows);
  - per site state: best-leaf purity (leaves >= MINN) and fragmentation F80 = the fewest leaves holding 80 % of the state's particles;
  - repeated splits: per finished node, eta2 of the site indicator over the node's refined states and over its decided groups; a node
    SPLITS the site when the decided groups' eta2 >= SPLIT_ETA (the number of such nodes beyond the first = repeated splits);
  - combination coverage (discrete sets): combinations with >= 5000 particles that own a leaf >= 70 % pure;
  - site-local map fidelity (--maps): per leaf map (mean of its unfiltered halves, lp LP, scaled to the truth on the molecule outside the
    sites) the site OCCUPANCY inside the site mask, occ = <x - t0, t1 - t0>_M / |t1 - t0|_M^2 (t0 / t1 = truth of the leaf's majority
    combination with the site in state 0 / 1), against the leaf's true state-1 share; and the masked correlation with the majority truth.
An unfinished tree is scored at its current state: finished nodes' leaves plus every pending node as one pseudo-leaf.
Continuous sites (simarmsigm1 angles) use the sign of the angle as the binary state for the discrete metrics.
usage: score_tree_sites.py --tree DIR --truth TRUTH [--maps DIR --masks A,B,C --minn 500 --split_eta 0.1 --lp 10 --json OUT]"""
import os, glob, json, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('--tree', required=True); ap.add_argument('--truth', required=True); ap.add_argument('--maps', default=''); ap.add_argument('--masks', default='')
ap.add_argument('--minn', type=int, default=500); ap.add_argument('--split_eta', type=float, default=0.1); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--json', default='')
ap.add_argument('--combo_min', type=int, default=5000); ap.add_argument('--pure', type=float, default=0.7)
ap.add_argument('--labels', default='', help='score this labelling (per master row class id, 0 = none) instead of the tree partition'); ap.add_argument('--leafmaps', default='', help='with --labels: lines "name even odd" per class id 1..K')
a = ap.parse_args(); R = a.tree
hdr = [l.split()[1:] for l in open(a.truth) if l.startswith('#')]; names = hdr[0][3:]; kinds = hdr[1][1:]
T = np.genfromtxt(a.truth, comments='#'); N = len(T); combo = T[:, 1].astype(int); junk = T[:, 2] > 0; X = T[:, 3:]
Xb = np.where(np.isfinite(X), (X > 0).astype(float), np.nan)   # binary state per site (discrete: as is; continuous: sign)
def nmi(x, y):
    _, xi = np.unique(x, return_inverse=True); _, yi = np.unique(y, return_inverse=True); J = np.zeros((xi.max() + 1, yi.max() + 1)); np.add.at(J, (xi, yi), 1); J /= J.sum()
    px, py = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(px, py)[nz])).sum(); hx, hy = -(px * np.log(px)).sum(), -(py * np.log(py)).sum(); return I / max(np.sqrt(hx * hy), 1e-12)
def eta2(v, g):
    m = np.isfinite(v); v, g = v[m], g[m]
    if len(v) < 10 or v.var() == 0: return np.nan
    _, gi = np.unique(g, return_inverse=True); cnt = np.bincount(gi); mu = np.bincount(gi, v) / cnt; return float((cnt * (mu - v.mean()) ** 2).sum() / (len(v) * v.var()))
def rows_of(lab_file, sid): return np.loadtxt(lab_file, dtype=int) == int(sid)
# ---- current partition: finished nodes' leaves + pending nodes
nodes = sorted(d[len(R) + 6:] for d in glob.glob(R + '/node_*') if os.path.isdir(d)); done = [u for u in nodes if os.path.exists('%s/node_%s/node.done' % (R, u))]
leaves = []   # (name, rows bool, map triple or None, junkflag, node)
if a.labels:
    fl = np.loadtxt(a.labels, dtype=int); LM = [l.split() for l in open(a.leafmaps) if l.strip() and not l.startswith('#')] if a.leafmaps else []
    for i in range(1, fl.max() + 1): leaves.append((LM[i - 1][0] if i <= len(LM) else 'class%d' % i, fl == i, (None,) + tuple(LM[i - 1][1:3]) if i <= len(LM) else None, False, 'labels'))
    state = 'labelling %s' % a.labels; done = []
elif os.path.exists(R + '/final_labels.txt') and os.path.exists(R + '/leaves.txt'):
    fl = np.loadtxt(R + '/final_labels.txt', dtype=int)
    for i, l in enumerate([l.split() for l in open(R + '/leaves.txt') if l.strip()], 1): leaves.append((l[0], fl == i, l[1:4], len(l) > 7 and l[7] == 'junk', l[0].split('_')[0]))
    state = 'finished'
else:
    for u in done:
        nl = '%s/node_%s/node_leaves.txt' % (R, u)
        for l in [l.split() for l in open(nl) if l.strip()] if os.path.exists(nl) else []: leaves.append((l[0], rows_of(l[5], l[6]), l[1:4], len(l) > 7 and l[7] == 'junk', u))
        ch = '%s/node_%s/children.txt' % (R, u)
        for l in [l.split() for l in open(ch) if l.strip()] if os.path.exists(ch) else []:
            if l[0] not in done: leaves.append(('pending_' + l[0], np.loadtxt(l[3], dtype=int) > 0, None, False, l[0]))
    if not done: leaves.append(('pending_r', np.ones(N, bool), None, False, 'r'))
    state = 'in progress (%d of %d node dirs finished)' % (len(done), len(nodes))
lab = np.zeros(N, int)
for i, L in enumerate(leaves, 1): lab[L[1]] = i
cov = lab > 0; ok = cov & ~junk; sizes = np.bincount(lab, minlength=len(leaves) + 1)
out = {'tree': R, 'state': state, 'leaves': len(leaves), 'pending': sum(1 for L in leaves if L[2] is None), 'sites': {}}
print('== %s: %s; %d leaves (%d pending pseudo-leaves, %d flagged junk), %d of %d rows covered' % (R, state, len(leaves), out['pending'], sum(L[3] for L in leaves), cov.sum(), N))
jl = np.array([False] + [L[3] for L in leaves]); tj = junk & cov
if junk.any(): print('   junk: %d true junk rows; in junk-flagged leaves %d (recall %.3f), junk-flagged rows %d of which true junk %.3f' % (junk.sum(), (tj & jl[lab]).sum(), (tj & jl[lab]).sum() / max(junk.sum(), 1), jl[lab].sum(), (tj & jl[lab]).sum() / max(jl[lab].sum(), 1)))
if (combo > 0).any() and kinds[0] == 'd': print('   NMI(combination; leaf) %.3f over %d non-junk rows' % (nmi(combo[ok], lab[ok]), ok.sum()))
big = [i for i in range(1, len(leaves) + 1) if sizes[i] >= a.minn]
for j, nm in enumerate(names):
    v, vb = X[:, j], Xb[:, j]; m = ok & np.isfinite(vb); rec = {'nmi': nmi(vb[m], lab[m]), 'eta2': eta2(v[m], lab[m]) if kinds[j] == 'c' else eta2(vb[m], lab[m])}
    line = '-- site %s: NMI(state; leaf) %.3f, eta2 %.3f' % (nm, rec['nmi'], rec['eta2'])
    for st in (1, 0):
        s_ = m & (vb == st); tot = s_.sum()
        if tot == 0: continue
        cnt = np.bincount(lab[s_], minlength=len(leaves) + 1); pur = [(cnt[i] / sizes[i], i) for i in big]; bp, bi = max(pur) if pur else (0., 0)
        o = np.sort(cnt[1:])[::-1]; f80 = int(np.searchsorted(np.cumsum(o), 0.8 * tot) + 1)
        rec['state%d' % st] = {'n': int(tot), 'share': float(tot / m.sum()), 'best_purity': float(bp), 'best_leaf': leaves[bi - 1][0] if bi else '', 'best_leaf_n': int(sizes[bi]) if bi else 0, 'F80': f80}
        line += '; state %d: %d (%.3f), best-leaf purity %.2f (%s, %d ptcls), F80 %d leaves' % (st, tot, tot / m.sum(), bp, leaves[bi - 1][0] if bi else '-', sizes[bi] if bi else 0, f80)
    out['sites'][nm] = rec; print(line)
# ---- per-node splits
if done: print('-- per finished node: eta2 of each site indicator over the refined states | over the decided groups (* = split, eta2 >= %.2f)' % a.split_eta)
splits = {nm: [] for nm in names}
for u in done:
    D = '%s/node_%s' % (R, u); lf = D + '/refined_s.txt' if os.path.exists(D + '/refined_s.txt') else D + '/refined_f.txt'
    if not os.path.exists(lf) or not os.path.exists(D + '/groups.txt'): continue
    st = np.loadtxt(lf, dtype=int); grp = np.zeros(N, int); gl = [l.split() for l in open(D + '/groups.txt') if l.strip()]
    for g, l in enumerate(gl, 1): grp[np.isin(st, [int(x) for x in l[2].split('+')]) & (st > 0)] = g
    mm = (st > 0) & ~junk; ng = len(set(grp[mm & (grp > 0)])); cells = []
    for j, nm in enumerate(names):
        vb = Xb[:, j]; m2 = mm & np.isfinite(vb); e_s = eta2(vb[m2], st[m2]); e_g = eta2(vb[m2 & (grp > 0)], grp[m2 & (grp > 0)]) if ng >= 2 else 0.
        sp = ng >= 2 and e_g >= a.split_eta
        if sp: splits[nm].append(u)
        cells.append('%s %.3f|%.3f%s' % (nm, e_s, e_g, '*' if sp else ''))
    print('   node %-6s %7d ptcls, %d states, %d groups (%s): %s' % (u, mm.sum(), len(set(st[st > 0])), ng, ' '.join(l[0][0] for l in gl), '  '.join(cells)))
for nm in names:
    out['sites'][nm]['split_nodes'] = splits[nm]
    if done: print('   site %s split at %d node(s) %s -> repeated splits %d' % (nm, len(splits[nm]), ' '.join(splits[nm]), max(0, len(splits[nm]) - 1)))
# ---- combination coverage
if kinds[0] == 'd':
    cc = [c for c in np.unique(combo[ok]) if (combo[ok] == c).sum() >= a.combo_min]; covd = []
    for c in np.unique(combo[ok]):
        cnt = np.bincount(lab[ok & (combo == c)], minlength=len(leaves) + 1); pur = max([(cnt[i] / sizes[i], i) for i in big] or [(0., 0)])
        print('   combination %d (%s): %6d ptcls, best leaf purity %.2f (%s)' % (c, ' '.join('%s=%d' % (nm[0], X[np.where(combo == c)[0][0], j]) for j, nm in enumerate(names)), (combo[ok] == c).sum(), pur[0], leaves[pur[1] - 1][0] if pur[1] else '-'))
        if c in cc and pur[0] >= a.pure: covd.append(int(c))
    out['coverage'] = {'big_combos': [int(c) for c in cc], 'covered': covd}; print('-- combination coverage: %d of %d combinations with >= %d particles own a leaf >= %.0f %% pure %s' % (len(covd), len(cc), a.combo_min, 100 * a.pure, covd))
# ---- site-local map fidelity
if a.maps and a.masks:
    def rd(p):
        h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
    def fcrop(v, nb):
        n = v.shape[0]
        if n == nb: return v
        F = np.fft.fftshift(np.fft.fftn(v)); lo = n // 2 - nb // 2; return np.real(np.fft.ifftn(np.fft.ifftshift(F[lo:lo + nb, lo:lo + nb, lo:lo + nb]))) * (nb / n) ** 3
    withmap = [i for i in big if leaves[i - 1][2] is not None]
    if withmap:
        nb = int(np.fromfile(leaves[withmap[0] - 1][2][1], np.int32, 1)[0]); smpd = float(np.fromfile(leaves[withmap[0] - 1][2][1], np.float32, 13)[10]) / nb
        f = np.fft.fftfreq(nb, d=smpd); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(nb, d=smpd), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
        k1, k2_ = 1. / (a.lp + 1.), 1. / (a.lp - 0.5); filt = np.where(k <= k1, 1., np.where(k >= k2_, 0., 0.5 + 0.5 * np.cos(np.pi * (k - k1) / (k2_ - k1))))
        lpf = lambda v: np.fft.irfftn(np.fft.rfftn(v) * filt, s=v.shape)
        combos_ = sorted(set(int(c) for c in np.unique(combo[combo > 0]))); TM = {c: lpf(fcrop(rd('%s/conf%02d_sites.mrc' % (a.maps, c)), nb)) for c in combos_}
        MS = [np.clip(fcrop(rd(p), nb), 0, 1) for p in a.masks.split(',')]; allm = np.clip(sum(MS), 0, 1)
        g = np.arange(nb) - nb // 2; sph = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) <= nb * 0.45
        cmap = {int(c): X[np.where(combo == c)[0][0]].astype(int) for c in combos_}
        def combo_of(st): return [c for c in combos_ if (cmap[c] == st).all()][0]
        res = {nm: [] for nm in names}
        for i in withmap:
            L = leaves[i - 1]; x = lpf(0.5 * (rd(L[2][1]) + rd(L[2][2]))); rows = (lab == i) & ok
            if rows.sum() < a.minn: continue
            cm = np.bincount(combo[rows]).argmax(); t = TM[cm]; w = sph & (allm < 0.05) & (np.abs(t) > 0.1 * np.abs(t).max())
            sc = (x[w] @ t[w]) / max(x[w] @ x[w], 1e-30); xs = x * sc   # scale on the molecule outside the sites
            for j, nm in enumerate(names):
                st = cmap[cm].copy(); st[j] = 0; t0 = TM[combo_of(st)]; st[j] = 1; t1 = TM[combo_of(st)]; M_ = MS[j]; d = t1 - t0
                occ = float(((xs - t0) * d * M_).sum() / max((d * d * M_).sum(), 1e-30)); fr = float(np.nanmean(Xb[rows, j]))
                mc = float(np.corrcoef((xs * M_)[M_ > 0.05], (t * M_)[M_ > 0.05])[0, 1]); res[nm].append((L[0], int(rows.sum()), fr, occ, mc))
        for nm in names:
            r_ = np.array([(n_, fr, oc, mc) for _, n_, fr, oc, mc in res[nm]])
            if len(r_) < 2: continue
            wcor = np.corrcoef(r_[:, 1], r_[:, 2])[0, 1] if r_[:, 1].std() > 0 else np.nan; rmse = np.sqrt(np.average((r_[:, 2] - r_[:, 1]) ** 2, weights=r_[:, 0]))
            top = sorted(res[nm], key=lambda t: -t[2])[:3]
            out['sites'][nm]['map'] = {'occ_vs_share_corr': float(wcor), 'occ_rmse': float(rmse), 'masked_corr_median': float(np.median(r_[:, 3])), 'top': [(t[0], t[1], t[2], t[3]) for t in top]}
            print('-- site %s maps (%d leaves): occupancy vs true share corr %.3f, rmse %.3f; masked corr with the majority truth median %.3f; leaves richest in state 1: %s' % (nm, len(r_), wcor, rmse, np.median(r_[:, 3]), '; '.join('%s %d ptcls share %.2f occ %.2f' % t[:4] for t in top)))
if a.json: json.dump(out, open(a.json, 'w'), indent=1, default=float)

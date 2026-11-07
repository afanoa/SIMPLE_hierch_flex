#!/usr/bin/env python3
"""Post-hoc factorization of a hetero_tree's leaf maps into independent difference components (truth only scores).
usage: mm_leaf_factors.py TREE_DIR [lp=12] [smpd=160x1.464/box] [synth=yes] [truth=..] [factors=..] [masks=DIR] [out=DIR]
1 every leaf's even / odd unfiltered half-maps (hetero_tree_leaves.txt columns 5-6), Fourier low-passed to lp A, inside an envelope made
  from the size-weighted mean map (data-derived: the smoothed mean above 0.15 of its max, dilated);
2 the cross-half reproducible rank r of the size-weighted, centred leaf maps: eigenvalues of the symmetrised even x odd Gram above
  3 x |its most negative eigenvalue| (the noise spread);
3 spatial ICA over the envelope voxels of the rank-r subspace -> r difference components; each leaf's code = its LS coordinate on them;
4 factor-state groups without re-assignment: per component, the leaves 2-means-split by code (size-weighted) -> two unions of leaves;
5 dependence: per component pair, the size-weighted correlation of leaf codes (independent motions ~ 0);
6 truth (optional): each component's overlap with the truth domain masks, each leaf code's correlation with the leaf's truth factor
  fractions, and the pooled groups' factor purity and map contrast vs the best leaf pair. Writes OUT/component_k.mrc and codes.txt."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
from sklearn.decomposition import FastICA
args = dict(a.split('=', 1) if '=' in a else ('tree', a) for a in sys.argv[1:]); T = args['tree']; lp = float(args.get('lp', 12))
out = args.get('out', T.rstrip('/') + '/leaf_factors'); os.makedirs(out, exist_ok=True)
L = [l.split() for l in open(T + '/hetero_tree_leaves.txt') if l.strip()]
names = [l[0] for l in L]; n = np.array([float(l[2]) for l in L]); fe = [l[5] for l in L]; fo = [l[6] for l in L]
def loc(p): return p if os.path.exists(p) else os.path.join(T, os.path.relpath(p, p[:p.index('1_hetero_tree') + len('1_hetero_tree')]) if '1_hetero_tree' in p else p)
E = np.array([read_mrc(loc(f)).astype(np.float32) for f in fe]); O = np.array([read_mrc(loc(f)).astype(np.float32) for f in fo]); B = E.shape[1]
smpd = float(args.get('smpd', 160 * 1.464 / B))   # the leaf maps sit at the tree's crop box (88 for box 160 at 1.464 A)
f = np.fft.fftfreq(B, smpd); fr = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2); filt = (fr <= 1. / lp).astype(np.float32)
def lowp(v): return np.real(np.fft.ifftn(np.fft.fftn(v) * filt)).astype(np.float32)
E = np.array([lowp(v) for v in E]); O = np.array([lowp(v) for v in O])
mean = ((E + O) * n[:, None, None, None]).sum(0) / (2 * n.sum()); sm = gaussian_filter(mean, 2.); env = binary_dilation(sm > 0.15 * sm.max(), iterations=4)
print('%d leaves (%s rows), lp %.0f A, envelope %d voxels' % (len(names), int(n.sum()), lp, env.sum()))
w = n / n.sum(); Xe = E[:, env]; Xo = O[:, env]; Xe = Xe - (w[:, None] * Xe).sum(0); Xo = Xo - (w[:, None] * Xo).sum(0)
sw = np.sqrt(n)[:, None]; G = (sw * Xe) @ (sw * Xo).T; ev, U = np.linalg.eigh(0.5 * (G + G.T)); ev, U = ev[::-1], U[:, ::-1]
thr = 3 * abs(min(ev.min(), 0)); r = int((ev > thr).sum()) if 'rank' not in args else int(args['rank'])   # selection imprint keeps every cross-half eigenvalue positive on real trees: rank= overrides
print('cross-half Gram eigenvalues (relative to the largest): %s; noise threshold %.3g -> rank %d' % (' '.join('%.3f' % (e / ev[0]) for e in ev[:8]), thr / ev[0], r))
X = sw * 0.5 * (Xe + Xo); Vt = np.linalg.svd(X, full_matrices=False)[2][:max(r, 1)]
if r >= 2: ica = FastICA(n_components=r, whiten='unit-variance', random_state=0, max_iter=4000); D = ica.fit_transform(Vt.T).T
else: D = Vt[:1]
D /= np.linalg.norm(D, axis=1, keepdims=True); Xm = 0.5 * (Xe + Xo); C = np.linalg.lstsq(D.T, Xm.T, rcond=None)[0].T   # leaves x r codes
for k in range(len(D)):
    sk = np.sign(D[k][np.argmax(np.abs(D[k]))]); D[k] *= sk; C[:, k] *= sk; v = np.zeros_like(mean); v[env] = D[k]; write_mrc('%s/component_%d.mrc' % (out, k + 1), v, smpd)
groups = []
for k in range(len(D)):
    t = C[:, k]; c = np.percentile(t, [25, 75])
    for _ in range(50):
        g = np.abs(t - c[1]) < np.abs(t - c[0])
        if g.all() or (~g).all(): break
        c = np.array([np.average(t[~g], weights=n[~g]), np.average(t[g], weights=n[g])])
    groups.append(g)
with open(out + '/codes.txt', 'w') as fw:
    fw.write('# leaf rows ' + ' '.join('c%d' % (k + 1) for k in range(len(D))) + ' ' + ' '.join('g%d' % (k + 1) for k in range(len(D))) + '\n')
    for i in range(len(names)): fw.write('%s %d %s %s\n' % (names[i], n[i], ' '.join('%.4f' % x for x in C[i]), ' '.join(str(int(g[i]) + 1) for g in groups)))
if len(D) > 1:
    cw = np.cov(C.T, aweights=n); cr = cw / np.sqrt(np.outer(np.diag(cw), np.diag(cw)))
    print('leaf-code correlation between components (size-weighted): max |r| %.2f\n%s' % (np.abs(cr - np.eye(len(D))).max(), np.array2string(cr, precision=2)))
print('component: leaves in group 2 (rows)')
for k in range(len(D)): print('  c%d: %s (%d rows)' % (k + 1, ' '.join(names[i] for i in np.where(groups[k])[0]), n[groups[k]].sum()))
# truth
if 'factors' in args:
    md = args.get('masks', os.path.dirname(args['factors'])); M = []
    k = 1
    from scipy.ndimage import zoom
    while os.path.exists('%s/domain%d_mask.mrc' % (md, k)):
        mk = read_mrc('%s/domain%d_mask.mrc' % (md, k)).astype(np.float32); M.append((zoom(mk, B / mk.shape[0], order=1) if mk.shape[0] != B else mk) > 0.5); k += 1
    print('truth: component energy share inside each domain mask (rows: components)')
    for j in range(len(D)):
        v = np.zeros_like(mean); v[env] = D[j]; e = v ** 2; print('  c%d: %s  outside %.2f' % (j + 1, ' '.join('d%d %.2f' % (q + 1, e[m].sum() / e.sum()) for q, m in enumerate(M)), 1 - sum(e[m].sum() for m in M) / e.sum()))
    t = np.loadtxt(args['truth'], ndmin=2); rows = t[:, 0].astype(int); combo = t[:, 1].astype(int); junk = t[:, 2] > 0
    F = np.loadtxt(args['factors'], ndmin=2); nf = F.shape[1] - 3; fmap = {int(x[0]): x[1:1 + nf].astype(int) for x in F}
    lab = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1)[rows - 1]; gen = ~junk
    fs = np.zeros((len(rows), nf), int)
    for i in np.where(gen)[0]: fs[i] = fmap[combo[i]]
    frac = np.array([[(fs[(lab == i + 1) & gen, q] == 2).mean() if ((lab == i + 1) & gen).any() else np.nan for q in range(nf)] for i in range(len(names))])
    ok = ~np.isnan(frac).any(1)
    print('truth: corr(leaf code, leaf fraction in state 2) (rows: components, cols: factors), leaves with genuine rows')
    for j in range(len(D)): print('  c%d: %s' % (j + 1, ' '.join('f%d %+.2f' % (q + 1, np.corrcoef(C[ok, j], frac[ok, q])[0, 1]) for q in range(nf))))
    print('truth: pooled factor-state groups (union of leaves by code) vs the best leaf pair, genuine rows')
    for q in range(nf):
        st = fs[:, q] == 2; best = None
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                a_, b_ = gen & (lab == i + 1), gen & (lab == j + 1)
                if a_.sum() and b_.sum():
                    c = abs(st[a_].mean() - st[b_].mean()) * np.sqrt(a_.sum() * b_.sum() / (a_.sum() + b_.sum()))
                    if best is None or c > best[0]: best = (c, names[i], names[j], st[a_].mean(), st[b_].mean(), a_.sum(), b_.sum())
        res = []
        for j in range(len(D)):
            g2 = gen & np.isin(lab, np.where(groups[j])[0] + 1); g1 = gen & (lab > 0) & ~np.isin(lab, np.where(groups[j])[0] + 1)
            if g1.sum() and g2.sum(): res.append((abs(st[g2].mean() - st[g1].mean()) * np.sqrt(g1.sum() * g2.sum() / (g1.sum() + g2.sum())), j, st[g1].mean(), st[g2].mean(), g1.sum(), g2.sum()))
        bj = max(res) if res else None
        print('  f%d: best leaf pair %s/%s frac %.2f/%.2f n %d/%d contrast %.0f | pooled c%d groups frac %.2f/%.2f n %d/%d contrast %.0f' % ((q + 1,) + best[1:3] + best[3:5] + best[5:7] + (best[0],) + ((bj[1] + 1, bj[2], bj[3], bj[4], bj[5], bj[0]) if bj else (0, 0, 0, 0, 0, 0))))
# 7 split directions: every split's child-difference map (children = size-weighted means of their descendant leaves), its loadings on the
# components (LS) and the cross-half cosine between splits (even of one x odd of the other: noise-unbiased); same motion in two
# branches = a shared sub-split (the DAG / lattice question)
def node_of(nm): return nm.split('_')[0]
nodes = sorted({node_of(x)[:d] for x in names for d in range(1, len(node_of(x)) + 1)}, key=lambda x: (len(x), x)); dirs = []
for nd in nodes:
    ch = {}
    for i, nm in enumerate(names):
        p = node_of(nm)
        if not p.startswith(nd): continue
        key = p[len(nd)] if len(p) > len(nd) else (nm.split('_', 1)[1] if p == nd and '_' in nm else None)
        if key is not None: ch.setdefault(key, []).append(i)
    if len(ch) < 2: continue
    cm = {k: ((n[v][:, None] * Xe[v]).sum(0) / n[v].sum(), (n[v][:, None] * Xo[v]).sum(0) / n[v].sum(), n[v].sum()) for k, v in sorted(ch.items())}
    ks = sorted(cm); tot = sum(cm[k][2] for k in ks)
    if len(ks) == 2: pairs = [(ks[0], ks[1])]
    else: pairs = [(k, None) for k in ks]   # K-way: each child against the node mean
    for k1, k2 in pairs:
        if k2 is None:
            me = sum(cm[k][0] * cm[k][2] for k in ks) / tot; mo = sum(cm[k][1] * cm[k][2] for k in ks) / tot; de, do = cm[k1][0] - me, cm[k1][1] - mo; lab_ = '%s:%s-mean' % (nd, k1)
        else: de, do = cm[k1][0] - cm[k2][0], cm[k1][1] - cm[k2][1]; lab_ = '%s:%s-%s' % (nd, k1, k2)
        dirs.append((lab_, de, do))
print('split directions: loadings on the components (LS, half-average) and cross-half norm')
for lab_, de, do in dirs:
    co = np.linalg.lstsq(D.T, 0.5 * (de + do), rcond=None)[0]; nn = np.sqrt(max(de @ do, 0)); res = 0.5 * (de + do) - co @ D
    print('  %-14s |d| %.3g  %s  residual share %.2f' % (lab_, nn, ' '.join('c%d %+.3g' % (k + 1, x) for k, x in enumerate(co)), (res @ res) / max(0.25 * (de + do) @ (de + do), 1e-30)))
if len(dirs) > 1:
    print('split directions: cross-half cosine (even_a . odd_b / norms)')
    nrm = [np.sqrt(max(de @ do, 1e-30)) for _, de, do in dirs]
    for i, (la, dea, doa) in enumerate(dirs):
        print('  %-14s %s' % (la, ' '.join('%+.2f' % (0.5 * (dea @ dob + deb @ doa) / (nrm[i] * nrm[j])) for j, (lb, deb, dob) in enumerate(dirs))))
# 8 synth=yes: 2^r full-resolution combination maps for a factorial re-assignment (one stage-B-like refine3D pass), all from the data:
# D_full = LS regression of the centred full-resolution leaf maps (even + odd) on the leaf codes; per component the two levels are the
# 5th / 95th percentile of the particles' leaf codes (the purest leaves, not the diluted group means); initial state per project row from
# its leaf's group bits (0 for rows without a leaf). Writes OUT/comb_NN.mrc, OUT/comb_table.txt, OUT/comb_states.txt
if args.get('synth', 'no') == 'yes':
    r = len(D); Ef = np.array([read_mrc(loc(x)).astype(np.float32) for x in fe]); Of = np.array([read_mrc(loc(x)).astype(np.float32) for x in fo])
    Mf = 0.5 * (Ef + Of).reshape(len(names), -1); del Ef, Of
    mu_f = (w[:, None] * Mf).sum(0); Cc = C - (w[:, None] * C).sum(0)
    Dfull = np.linalg.lstsq(np.sqrt(n)[:, None] * Cc, np.sqrt(n)[:, None] * (Mf - mu_f), rcond=None)[0]   # r x V
    labr = np.loadtxt(T + '/hetero_tree_labels.txt', dtype=int, ndmin=1)
    pc = np.concatenate([np.full(int(n[i]), i) for i in range(len(names))])
    lv = np.array([np.percentile(Cc[pc, k], [5, 95]) for k in range(r)])
    with open(out + '/comb_table.txt', 'w') as ft:
        ft.write('# state ' + ' '.join('c%d_level' % (k + 1) for k in range(r)) + '  (1 = low code, 2 = high code; levels %s)\n' % ' '.join('%.4g/%.4g' % tuple(x) for x in lv))
        for st in range(2 ** r):
            bits = [(st >> k) & 1 for k in range(r)]; v = mu_f + sum(lv[k, bits[k]] * Dfull[k] for k in range(r))
            write_mrc('%s/comb_%02d.mrc' % (out, st + 1), v.reshape(E.shape[1:]).astype(np.float32), smpd); ft.write('%d %s\n' % (st + 1, ' '.join(str(b + 1) for b in bits)))
    gbit = np.array([sum(int(groups[k][i]) << k for k in range(r)) + 1 for i in range(len(names))])
    st0 = np.where(labr > 0, gbit[np.maximum(labr, 1) - 1], 0); rs = np.random.default_rng(5)
    for s_ in range(1, 2 ** r + 1):   # refine3D's multi-state startup needs every state populated: 50 random rows of the largest state
        if (st0 == s_).sum() < 50:
            big = np.bincount(st0[st0 > 0]).argmax(); st0[rs.choice(np.where(st0 == big)[0], 50, replace=False)] = s_
    np.savetxt(out + '/comb_states.txt', st0, fmt='%d')
    print('synth: %d combination maps, levels per component %s; initial states %s' % (2 ** r, ' '.join('%.3g/%.3g' % tuple(x) for x in lv), ' '.join('%d:%d' % (s_, (st0 == s_).sum()) for s_ in range(1, 2 ** r + 1))))

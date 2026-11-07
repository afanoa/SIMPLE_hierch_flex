#!/usr/bin/env python3
"""mest_prior.py TREE OUT [--truth T --csv C --maps M] [--closed] [--lp 10]
Final-map estimators of a finished hetero_tree's leaf units from the tree's own raw leaf halves
(TREE/final_maps/raw/halves/unit_NN_{even,odd}.mrc: masked, unfiltered), re-implementing tree_wiener's per-shell cross-half
covariance Wiener (two parts: inside / outside the focus mask where the units differ) and adding the tree-prior variants:
  raw     the plain leaf maps                          own     per-unit Wiener d / (d + nv)
  emp2    PSD empirical covariance (tree_wiener's delivered estimator); emp1 = the same without the focus-mask split
  clu2    tree covariance (NNLS) on a complete-linkage clustering of the units (tree_wiener's tree2)
  topo2   tree covariance (NNLS) on the hetero_tree's own topology (root, nodes, leaves)
  depth2  depth-tied tree covariance: cov(p, q) = sum_{d <= depth(lca)} theta_d, one variance per depth per shell
  anc2    sequential hierarchical MAP: child = parent estimate + Wiener-shrunk innovation, innovation variance per depth
Writes OUT/<est>/unit_NN.mrc (merged, low-passed) and OUT/<est>/halves/unit_NN_{even,odd}.mrc, OUT/summary.txt (covariance
R^2 per model, per-unit half-map FSC 0.143), the qtoy truth readout (mest_score.py; --truth/--csv/--maps) and the RNase P
closed-arm / open-arm occupancies (--closed)."""
import argparse, os, sys, glob, subprocess
import numpy as np
from scipy.optimize import nnls
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

H = '/mnt/beegfs/elmlund/afan/hetero_tree_20261028'
SC = H + '/fixtures/scripts'; MEST_SCORE = H + '/softhard/scripts/mest_score.py'; R160 = '/mnt/beegfs/elmlund/afan/routes/toy55/maps160'
REF = '/mnt/beegfs/elmlund/afan/hetero_tree_regression_20261026/runs/real_root_split_contrast_m8c52_s1/closed_reference_rows.txt'
HPFIT, LPFIT, R2MIN, LPMASK, MTHR, MSIG, MDIL, EDGE = 40., 8., 0.7, 20., 0.15, 2., 3., 6.
ESTS = ['raw', 'own', 'emp1', 'emp2', 'clu2', 'topo2', 'depth2', 'anc2']


def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); smpd = float(np.fromfile(p, np.float32, 13)[10]) / n
    return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64), smpd


def wr(p, x, smpd):
    h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd
    hf[13:16] = 90.; h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.astype(np.float32).tobytes())


def log(*a): print(*a, flush=True)


# ---------------------------------------------------------------- tree designs (tree_wiener's `designs`, generalised)

def design_from_parents(par, root, U):
    """par: id -> parent id for every node but the root; units are ids 0..U-1. Columns: root, then every other node; a pair's row
    has ones at the root and at the nodes shared by both units' paths (a unit's own node counts when p == q)."""
    paths = []
    for u in range(U):
        p, path = u, []
        while p != root:
            path.append(p); p = par[p]
        paths.append(path)
    ids = sorted(set(par) - {root})
    col = {x: 1 + j for j, x in enumerate(ids)}
    A = np.zeros((U * (U + 1) // 2, 1 + len(ids)))
    r = 0
    for p in range(U):
        for q in range(p, U):
            A[r, 0] = 1.
            for x in paths[p]:
                if x in paths[q]: A[r, col[x]] = 1.
            r += 1
    return A


def design_depth(depth_lca, U, D):
    """cov(p, q) = sum_{d <= depth(lca(p, q))} theta_d; depth_lca[p][q] with depth_lca[p][p] = the leaf's own depth."""
    A = np.zeros((U * (U + 1) // 2, D + 1)); r = 0
    for p in range(U):
        for q in range(p, U):
            A[r, :depth_lca[p][q] + 1] = 1.; r += 1
    return A


def pairs_to_mat(y, U):
    M = np.zeros((U, U)); r = 0
    for p in range(U):
        for q in range(p, U):
            M[p, q] = M[q, p] = y[r]; r += 1
    return M


def treefit(A, ce, U):
    y = np.array([ce[p, q] for p in range(U) for q in range(p, U)])
    th, _ = nnls(A, y); yf = A @ th
    off = np.array([p != q for p in range(U) for q in range(p, U)])
    ym = y[off].mean(); sres = ((y[off] - yf[off]) ** 2).sum(); stot = ((y[off] - ym) ** 2).sum()
    return pairs_to_mat(yf, U), 1. - sres / max(stot, 1e-300), th


def psd(ce):
    lam, V = np.linalg.eigh(ce); return (V * np.maximum(lam, 0.)[None]) @ V.T


def wiener_mats(cf, nv):
    W1 = cf @ np.linalg.inv(cf + np.diag(nv)); W2 = cf @ np.linalg.inv(cf + 2. * np.diag(nv))
    return W1, W2


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tree'); ap.add_argument('out'); ap.add_argument('--truth', default=''); ap.add_argument('--csv', default='')
    ap.add_argument('--maps', default=''); ap.add_argument('--closed', action='store_true'); ap.add_argument('--lp', type=float, default=10.)
    ap.add_argument('--ests', default=','.join(ESTS)); ap.add_argument('--closed_top', type=int, default=0, help='summary over the N leaves richest in reference closed rows')
    ap.add_argument('--halves_only', action='store_true', help='write no merged maps (large trees)')
    a = ap.parse_args(); T, out = a.tree, a.out; os.makedirs(out, exist_ok=True)
    ESTS[:] = [k for k in a.ests.split(',') if k]
    units = [l.split() for l in open(os.path.join(T, 'final_maps', 'units_index.txt')) if not l.startswith('#') and l.strip()]
    U = len(units); names = [int(u[0]) for u in units]; leaf = [u[1] for u in units]; counts = np.array([float(u[2]) for u in units]); ljunk = [u[3] == 'junk' for u in units]
    halves = os.path.join(T, 'final_maps', 'raw', 'halves')
    E, O = [], []
    for nm in names:
        e, smpd = rd(os.path.join(halves, 'unit_%02d_even.mrc' % nm)); o, _ = rd(os.path.join(halves, 'unit_%02d_odd.mrc' % nm)); E.append(e); O.append(o)
    n = E[0].shape[0]; nsh = n // 2; wl = counts / counts.sum()
    log('%s: %d units, box %d at %.3f A/px, counts %s' % (T, U, n, smpd, counts.astype(int).tolist()))
    # topology of the hetero_tree: leaf 'r4_s2' -> node 'r4' -> ... -> 'r'
    lnode = [l.split('_')[0] for l in leaf]; depth = [len(x) - 1 for x in lnode]
    # Fourier layout (rfftn): shells 1..nsh-1, Friedel weights 2 except the self-conjugate planes
    f = np.fft.fftfreq(n) * n; kz = np.arange(n // 2 + 1)
    kx, ky, kzz = np.meshgrid(f, f, kz, indexing='ij'); s3 = np.rint(np.sqrt(kx ** 2 + ky ** 2 + kzz ** 2)).astype(int)
    wg3 = np.where((kzz == 0) | (kzz == n // 2), 1., 2.)
    sel = (s3 >= 1) & (s3 < nsh); shell = s3[sel]; wg = wg3[sel]; nvox = int(sel.sum())
    fq = np.sqrt(kx ** 2 + ky ** 2 + kzz ** 2)[sel] / (n * smpd)
    k1, k2 = 1. / (a.lp + 1.), 1. / (a.lp - 0.5); filt = np.where(fq >= k2, 0., np.where(fq > k1, 0.5 + 0.5 * np.cos(np.pi * (fq - k1) / (k2 - k1)), 1.))
    order = np.argsort(shell, kind='stable'); shell_o = shell[order]; lo = np.searchsorted(shell_o, np.arange(nsh)); hi = np.searchsorted(shell_o, np.arange(nsh), 'right')
    sw = np.bincount(shell, wg, nsh)
    rsh = np.where(np.arange(nsh) > 0, n * smpd / np.maximum(np.arange(nsh), 1), np.inf)
    # focus mask: where the count-weighted cross-half product of the 20 A maps exceeds the product of the means
    def lpg(v, lp):
        F = np.fft.rfftn(v); f2 = (kx ** 2 + ky ** 2 + kzz ** 2) / float(n * n)
        return np.fft.irfftn(F * np.exp(-f2 / smpd ** 2 * lp ** 2 * np.log(2.)), s=v.shape)
    c = n // 2; zz, yy, xx = np.mgrid[:n, :n, :n]; rr = np.sqrt((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2)
    mskdiam = 0.   # the leaf nodes' plain reconstructions carry the tree's mskdiam (the root has none)
    for l in open(sorted(glob.glob(os.path.join(T, 'node_*', 'plain', 'hetero_step.sh')) + [os.path.join(T, 'node_r', 'fflex', 'hetero_step.sh')])[0]):
        if 'mskdiam=' in l: mskdiam = float(l.split('mskdiam=')[1].split()[0]); break
    r0 = mskdiam / 2. / smpd; msk = np.where(rr <= r0, 1., np.where(rr <= r0 + EDGE, 0.5 + 0.5 * np.cos(np.pi * (rr - r0) / EDGE), 0.)); inm = msk > 0.5
    eb = sum(wl[u] * lpg(E[u], LPMASK) for u in range(U)); ob = sum(wl[u] * lpg(O[u], LPMASK) for u in range(U))
    vv = sum(wl[u] * lpg(E[u], LPMASK) * lpg(O[u], LPMASK) for u in range(U)) - eb * ob; vv = gaussian_filter(vv, MSIG)
    core = (vv > MTHR * vv[inm].max()) & inm
    mv = np.clip(gaussian_filter((distance_transform_edt(~core) <= MDIL).astype(float), MSIG), 0., 1.) * msk if core.any() else np.zeros_like(msk)
    log('focus mask %.1f %% of the sphere (mskdiam %.1f A)' % (100. * mv.sum() / msk.sum(), mskdiam))
    coefs = lambda x: np.fft.rfftn(x)[sel].astype(np.complex64)
    ef = np.stack([coefs(E[u]) for u in range(U)], 1); of = np.stack([coefs(O[u]) for u in range(U)], 1)
    ea = np.stack([coefs(E[u] * mv) for u in range(U)], 1); oa = np.stack([coefs(O[u] * mv) for u in range(U)], 1)
    del E, O
    # clustering tree of the units (band 40-10 A correlation of the merged maps, complete linkage)
    band = (rsh[shell] <= 40.) & (rsh[shell] >= 10.)
    Xb = (0.5 * (ef[band] + of[band])).astype(complex) * np.sqrt(wg[band])[:, None]; Cm = np.real(Xb.conj().T @ Xb); dd = np.sqrt(np.diag(Cm)); corr = Cm / np.outer(dd, dd); del Xb
    Z = linkage(squareform(np.clip(1. - corr, 0., 2.), checks=False), 'complete')
    par = {}
    for i, (p, q, _, _) in enumerate(Z): par[int(p)] = U + i; par[int(q)] = U + i
    A_clu = design_from_parents(par, 2 * U - 2, U)
    # the hetero_tree's own topology: ids = units, then the nodes
    nodes = sorted({x[:k] for x in lnode for k in range(1, len(x) + 1)}, key=lambda x: (len(x), x)); nid = {x: U + j for j, x in enumerate(nodes)}
    par = {u: nid[lnode[u]] for u in range(U)}
    for x in nodes:
        if x != 'r': par[nid[x]] = nid[x[:-1]]
    A_topo = design_from_parents(par, nid['r'], U)
    D = max(depth) + 1
    def lca_depth(p, q):
        if p == q: return depth[p] + 1
        x, y = lnode[p], lnode[q]; k = 0
        while k < min(len(x), len(y)) and x[k] == y[k]: k += 1
        return k - 1
    A_depth = design_depth([[lca_depth(p, q) for q in range(U)] for p in range(U)], U, D)
    # pooled-map rows for anc2: node x -> count-weighted mean of its leaves (the leaf's own row is a unit vector)
    pool = {x: np.array([counts[u] if lnode[u].startswith(x) else 0. for u in range(U)]) for x in nodes}
    for x in nodes: pool[x] /= pool[x].sum()
    chain = {u: [lnode[u][:k] for k in range(1, len(lnode[u]) + 1)] for u in range(U)}   # root ... node of the leaf

    def cov(e, o, w):
        """cross-half covariance of the units over one shell's voxels (w = Friedel weights of those voxels) and per-unit noise."""
        g = np.real((e * w[:, None]).conj().T @ o)
        nv = np.maximum((w[:, None] * np.abs(e - o) ** 2).sum(0) / (4. * w.sum()), 1e-300)
        return 0.5 * (g + g.T) / w.sum(), nv

    est = {k: [np.zeros((nvox, U), np.complex64) for _ in range(3)] for k in ESTS}
    r2 = {k: np.zeros(nsh) for k in ('clu2', 'topo2', 'depth2', 'emp2', 'anc2')}

    def anc_mats(e, o, w, ce, nv_unit):
        """sequential hierarchical MAP as a U x U matrix per shell: innovation covariance per depth from the pooled halves."""
        inn_tau, inn_nv = {}, {}
        xh = lambda pe, po: (w * (pe.real * po.real + pe.imag * po.imag)).sum() / w.sum()
        nh = lambda pe, po: max((w * np.abs(pe - po) ** 2).sum() / (4. * w.sum()), 1e-300)
        for x in nodes:
            if x == 'r': continue
            pe = e @ pool[x] - e @ pool[x[:-1]]; po = o @ pool[x] - o @ pool[x[:-1]]
            inn_tau.setdefault(len(x) - 1, []).append((pool[x] @ counts, xh(pe, po))); inn_nv[x] = nh(pe, po)
        for u in range(U):   # the leaf below its node (zero when the node has one leaf)
            x = lnode[u]; pe = e[:, u] - e @ pool[x]; po = o[:, u] - o @ pool[x]
            inn_tau.setdefault(len(x), []).append((counts[u], xh(pe, po))); inn_nv[('leaf', u)] = nh(pe, po)
        tau = {d: max(sum(c * v for c, v in L) / max(sum(c for c, v in L), 1e-300), 0.) for d, L in inn_tau.items()}
        # root: its own pooled map shrunk by its cross-half power
        pr = pool['r']; dr = max(float(pr @ ce @ pr), 0.); nr = float((pr ** 2) @ nv_unit)
        W1 = np.zeros((U, U)); W2 = np.zeros((U, U))
        for u in range(U):
            row1 = dr / (dr + nr) * pr; row2 = dr / (dr + 2. * nr) * pr
            for x in chain[u][1:]:
                dlt = pool[x] - pool[x[:-1]]; d = len(x) - 1
                row1 = row1 + tau[d] / (tau[d] + inn_nv[x]) * dlt; row2 = row2 + tau[d] / (tau[d] + 2. * inn_nv[x]) * dlt
            dlt = np.eye(U)[u] - pool[lnode[u]]; d = len(lnode[u])
            row1 = row1 + tau[d] / (tau[d] + inn_nv[('leaf', u)]) * dlt; row2 = row2 + tau[d] / (tau[d] + 2. * inn_nv[('leaf', u)]) * dlt
            W1[u] = row1; W2[u] = row2
        return W1, W2

    def shell_pass(s, e, o, part):
        """part 0: full maps (raw / own / emp1); part 1 / 2: inside / outside the focus mask (the two-part estimators)."""
        idx = order[lo[s]:hi[s]]; w = wg[idx]
        ee, oo = e[idx].astype(complex), o[idx].astype(complex); ce, nv = cov(ee, oo, w); x = 0.5 * (ee + oo)
        if part == 0:
            for aa in range(U):
                d = max(ce[aa, aa], 0.)
                if 'raw' in est: est['raw'][0][idx, aa] = x[:, aa]; est['raw'][1][idx, aa] = ee[:, aa]; est['raw'][2][idx, aa] = oo[:, aa]
                if 'own' in est: est['own'][0][idx, aa] = d / (d + nv[aa]) * x[:, aa]; est['own'][1][idx, aa] = d / (d + 2 * nv[aa]) * ee[:, aa]; est['own'][2][idx, aa] = d / (d + 2 * nv[aa]) * oo[:, aa]
            if 'emp1' in est:
                W1, W2 = wiener_mats(psd(ce), nv)
                est['emp1'][0][idx] = x @ W1.T; est['emp1'][1][idx] = ee @ W2.T; est['emp1'][2][idx] = oo @ W2.T
            return
        mats = {}
        for k, A in (('clu2', A_clu), ('topo2', A_topo), ('depth2', A_depth)):
            if k in est: mats[k], r2[k][s], _ = treefit(A, ce, U)
        if 'emp2' in est: mats['emp2'] = psd(ce)
        for k, cf in mats.items():
            W1, W2 = wiener_mats(cf, nv)
            est[k][0][idx] += x @ W1.T; est[k][1][idx] += ee @ W2.T; est[k][2][idx] += oo @ W2.T
        if 'anc2' in est:
            W1, W2 = anc_mats(ee, oo, w, ce, nv)
            est['anc2'][0][idx] += x @ W1.T; est['anc2'][1][idx] += ee @ W2.T; est['anc2'][2][idx] += oo @ W2.T

    for s in range(1, nsh):
        shell_pass(s, ef, of, 0)
        shell_pass(s, ea, oa, 1)
        shell_pass(s, ef - ea, of - oa, 2)
    fit = (rsh >= LPFIT) & (rsh <= HPFIT) & (np.arange(nsh) > 0)
    # note: the two-part R^2 is that of the inside pass's covariance fit, as logged per shell by the last call
    R2K = [k for k in ('clu2', 'topo2', 'depth2') if k in est]
    log('covariance fit R^2 (off-diagonal, 40-8 A, the outside-mask part): ' + ' '.join('%s %.3f' % (k, r2[k][fit].mean()) for k in R2K))

    def fsc_res(e1, o1):
        num = np.bincount(shell, wg * (e1 * np.conj(o1)).real, nsh); de = np.bincount(shell, wg * np.abs(e1) ** 2, nsh); do = np.bincount(shell, wg * np.abs(o1) ** 2, nsh)
        c = num / np.maximum(np.sqrt(de * do), 1e-300)
        for s in range(1, nsh):
            if c[s] < 0.143: return n * smpd / max(s - 1, 1)
        return n * smpd / (nsh - 1)

    def put(cv, fn):
        F = np.zeros((n, n, n // 2 + 1), complex); F[sel] = cv; wr(fn, np.fft.irfftn(F, s=(n, n, n)), smpd)

    res = {}
    for k in ESTS:
        os.makedirs(os.path.join(out, k, 'halves'), exist_ok=True)
        for u in range(U):
            stem = ('junk_' if ljunk[u] else '') + 'unit_%02d' % names[u]
            if not a.halves_only: put(est[k][0][:, u] * filt, os.path.join(out, k, stem + '.mrc'))
            put(est[k][1][:, u], os.path.join(out, k, 'halves', stem + '_even.mrc')); put(est[k][2][:, u], os.path.join(out, k, 'halves', stem + '_odd.mrc'))
            res[(k, u)] = fsc_res(est[k][1][:, u], est[k][2][:, u])
    with open(os.path.join(out, 'summary.txt'), 'w') as fo:
        fo.write('# %s: %d units; half-map FSC 0.143 (A) per estimator and unit; covariance R^2 (40-8 A): %s\n' % (T, U, ' '.join('%s %.3f' % (k, r2[k][fit].mean()) for k in R2K)))
        fo.write('%-8s' % 'est' + ''.join('%9s' % leaf[u] for u in range(U)) + '   median\n')
        for k in ESTS:
            fo.write('%-8s' % k + ''.join('%9.2f' % res[(k, u)] for u in range(U)) + '   %6.2f\n' % np.median([res[(k, u)] for u in range(U)]))
    log(open(os.path.join(out, 'summary.txt')).read())
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    if a.truth:
        arms = ['%s=%s:halves/unit_{k:02d}_{eo}.mrc' % (k, os.path.join(out, k)) for k in ESTS]
        arms += ['TW_%s=%s:halves/unit_{k:02d}_{eo}.mrc' % (e, os.path.join(T, 'final_maps', e)) for e in ('raw', 'tree2', 'emp2') if os.path.isdir(os.path.join(T, 'final_maps', e, 'halves'))]
        os.makedirs(os.path.join(out, 'score'), exist_ok=True)
        subprocess.run([sys.executable, MEST_SCORE, T, a.truth, a.csv, a.maps, os.path.join(out, 'score')] + arms, check=False, env=env)
    if a.closed:
        sd = os.path.join(out, 'closed'); os.makedirs(sd, exist_ok=True); spec = []
        for k in ESTS + ['TW_raw', 'TW_tree2', 'TW_emp2']:
            d = os.path.join(out, k, 'halves') if not k.startswith('TW_') else os.path.join(T, 'final_maps', k[3:], 'halves')
            for u in range(U):
                e = os.path.join(d, 'unit_%02d_even.mrc' % names[u])
                if os.path.exists(e): spec.append('%s__%02d %s %s' % (k, names[u], e, e.replace('_even', '_odd')))
        open(os.path.join(sd, 'spec.txt'), 'w').write('\n'.join(spec) + '\n')
        subprocess.run([sys.executable, os.path.join(SC, 'make_view.py'), os.path.join(sd, 'views'), os.path.join(sd, 'spec.txt'), '--box', '160'], check=True, env=env, stdout=subprocess.DEVNULL)
        CL = sorted(glob.glob(os.path.join(R160, 'docked', 'conf18_*_m160.mrc')))[0]
        subprocess.run([sys.executable, os.path.join(SC, 'arm_readout.py'), os.path.join(R160, 'toy_consensus.mrc'), CL, os.path.join(sd, 'arm.txt')] +
                       [os.path.join(sd, 'views', s.split()[0] + '.mrc') for s in spec], check=True, env=env, stdout=subprocess.DEVNULL)
        arm = {os.path.basename(l.split()[0])[:-4]: (float(l.split()[1]), float(l.split()[2])) for l in open(os.path.join(sd, 'arm.txt')) if not l.startswith('#')}
        show = list(range(U))
        if a.closed_top:   # the leaves richest in the reference closed rows (the tree's labels index the leaves by unit name)
            lab = np.loadtxt(os.path.join(T, 'hetero_tree_labels.txt'), dtype=int); rr_ = np.loadtxt(REF, dtype=int)
            ref = np.zeros(len(lab) + 1, bool); ref[rr_[rr_ <= len(lab)]] = True
            share = [ref[np.where(lab == names[u])[0] + 1].mean() if (lab == names[u]).any() else 0. for u in range(U)]
            show = sorted(range(U), key=lambda u: -share[u])[:a.closed_top]
            log('closed-row share of the shown leaves: ' + ' '.join('%s %.3f' % (leaf[u], share[u]) for u in show))
        log('closed-arm / open-arm occupancy per unit (closed corner * = >= 0.95 / <= 0.6); corners = over all %d units' % U)
        log('%-8s %7s %7s' % ('est', 'corners', 'FSCmed') + ''.join('%14s' % leaf[u] for u in show))
        for k in ESTS + ['TW_raw', 'TW_tree2', 'TW_emp2']:
            cells = []
            for u in show:
                t = '%s__%02d' % (k, names[u])
                cells.append('%5.2f/%4.2f%s' % (*arm[t], '*' if arm[t][0] >= 0.95 and arm[t][1] <= 0.6 else ' ') if t in arm else '-')
            corners = sum(1 for u in range(U) if '%s__%02d' % (k, names[u]) in arm and arm['%s__%02d' % (k, names[u])][0] >= 0.95 and arm['%s__%02d' % (k, names[u])][1] <= 0.6)
            fmed = np.median([res[(k, u)] for u in range(U)]) if k in ESTS else float('nan')
            log('%-8s %7d %7.2f' % (k, corners, fmed) + ''.join('%14s' % x for x in cells))


if __name__ == '__main__':
    main()

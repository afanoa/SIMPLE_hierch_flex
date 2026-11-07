#!/usr/bin/env python3
"""mest_charts.py TREE OUT [--truth KIND:ARGS] [--nnull 100 --nproc 8 --kmax 8 --min_group 200 --zsrc txt|ml]
The hetero_tree as a hierarchical mixture of local PPCA charts (every node's fflex fit) with globally coordinated latent
coordinates, and a density-gap split criterion in each node's deconvolved latent density.
1. Coordination. A child's particles are embedded in the child's and the parent's chart; the child -> parent affine map is fitted
   by least squares on them, weighted by the parent coordinates' posterior precision (flex_pca_embedding.bin, carried into the
   delivered frame by the least-squares rotation between the raw and the delivered coordinates), and composed up to the root: a
   particle's global coordinate is its deepest chart's coordinate carried into the root frame. Readouts on the qtoy fixtures: the
   ordering of the arm conformers by s1 along the unsupervised first axis of the arm particles (Spearman, overall and per conformer
   mean), the closed and holo separation from the arm (cross-half LDA AUROC and d'), for the composed global coordinates, the root
   chart alone and (within node) the local charts.
2. Split criterion. For every candidate split of a node (split_children_pre.txt groups) and pair of children: the particles'
   projections on the pair's mean-difference direction in the node's chart with their measurement variances (zsrc=txt: the
   delivered MAP coordinates with flex's noise-scaled precisions; zsrc=ml: the MAP coordinates un-shrunk by the map prior, with the
   raw precisions), deconvolved in 1-D by extreme deconvolution (tied-variance mixture, K by cross-half held-out likelihood).
   Gap depth D = 1 - valley / lower peak of the deconvolved density between the children's modes, its p-value against a
   uniform-continuum null (same deconvolved mean and variance, same noise, same pipeline); the same D from flex's own deconvolved
   prior projected on the direction (D_flex); the children's d' (floored when the direction carries no latent signal).
Writes OUT/charts.txt and OUT/charts_summary.json."""
import argparse, json, os, sys, time
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.environ.get('MEST_LATENT', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent'))
import latent_states_v4 as LS

log = LS.log


# ---------------------------------------------------------------- readers

def read_embedding(fname):
    """raw block of flex_pca_embedding.bin: MAP z (n, k), posterior precision (n, k, k); trailing block: noise scale, deconvolved z."""
    with open(fname, 'rb') as f:
        magic = f.read(8).decode(); ver, n, k, box = np.frombuffer(f.read(16), np.int32); f.read(4)
        assert magic == 'SIMPLFXC' and ver == 5, '%s: %s v%d' % (fname, magic, ver)
        pinds = np.frombuffer(f.read(4 * n), np.int32).copy()
        z = np.frombuffer(f.read(8 * n * k), np.float64).reshape(k, n).T.copy()
        f.read(8 * k + 24 * n)
        prec = np.frombuffer(f.read(8 * k * k * n), np.float64).reshape(n, k, k).copy()
        f.read(8)
        out = dict(n=int(n), k=int(k), pinds=pinds, z=z, prec=prec, noise_scale=1., zdec=None, precdec=None)
        if f.read(8) == b'SIMPLFXD':
            out['noise_scale'] = float(np.frombuffer(f.read(8), np.float64)[0]); ld = np.frombuffer(f.read(4), np.int32)[0]
            if ld == 1:
                nb, kb = np.frombuffer(f.read(8), np.int32)
                out['zdec'] = np.frombuffer(f.read(8 * nb * kb), np.float64).reshape(kb, nb).T.copy()
                out['precdec'] = np.frombuffer(f.read(8 * kb * kb * nb), np.float64).reshape(nb, kb, kb).copy()
    return out


def read_coords(fname):
    d = np.loadtxt(fname, ndmin=2)
    return dict(pind=d[:, 0].astype(int), eo=d[:, 1].astype(int), label=d[:, 2].astype(int), z=d[:, 5:])


def read_deconv_prior(fname):
    hdr = open(fname).readline()
    ns = float(hdr.split('noise_scale=')[1].split()[0]); nc = int(hdr.split('ncomp=')[1].split()[0])
    d = np.loadtxt(fname, comments='#', ndmin=2)
    return dict(w=d[:, 1], mu=d[:, 2:2 + nc], sig=d[:, 2 + nc:2 + 2 * nc], noise_scale=ns)


def read_map_prior(fname):
    d = np.loadtxt(fname, comments='#', ndmin=2)
    return d[:, 2]                                                         # prior precision per component


def read_ints(f):
    return np.loadtxt(f, dtype=int, ndmin=1) if os.path.exists(f) and os.path.getsize(f) > 0 else None


def load_truth(spec, N):
    kind, _, rest = spec.partition(':')
    if kind == 'qtoy':
        hdr = open(rest).readline().lstrip('#').split(); t = np.loadtxt(rest, ndmin=2)
        assert (t[:, 0] == np.arange(1, len(t) + 1)).all() and len(t) >= N
        T = {c: t[:N, hdr.index(c)] for c in hdr if c in ('conformer', 'junk', 's1', 'holo', 'combo')}
        T['arm'] = (T['junk'] < 0.5) & (T['holo'] < 0.5) & (T['conformer'] != 18) & np.isfinite(T['s1'])
        return kind, T
    if kind == 'closed':
        fxr, clr = rest.split(':'); ref_rows = np.loadtxt(clr, dtype=int, ndmin=1)
        rm = np.loadtxt(fxr, dtype=int, ndmin=1) if os.path.exists(fxr) else None
        ix = rm[:N] if rm is not None and len(rm) >= N else np.arange(1, N + 1)
        return kind, dict(closedref=np.isin(ix, ref_rows).astype(float))
    return 'none', {}


# ---------------------------------------------------------------- statistics

def spearman(a, b):
    ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def auroc(score, pos):
    o = np.argsort(score, kind='stable'); r = np.empty(len(score)); r[o] = np.arange(1, len(score) + 1)
    npos = pos.sum(); nneg = len(pos) - npos
    return float((r[pos].sum() - npos * (npos + 1) / 2.) / max(npos * nneg, 1))


def lda_crosshalf(X, pos, eo):
    """direction fitted on one half (LDA with a ridge), scored on the other and vice versa: AUROC and d' of the pooled scores."""
    sc = np.zeros(len(X))
    for fs, es in ((eo == 0, eo == 1), (eo == 1, eo == 0)):
        a, b = X[fs & pos], X[fs & ~pos]
        if len(a) < 5 or len(b) < 5: return float('nan'), float('nan')
        S = 0.5 * (np.cov(a.T) + np.cov(b.T)) + 1e-3 * np.eye(X.shape[1]) * np.trace(np.cov(X[fs].T)) / X.shape[1]
        sc[es] = X[es] @ np.linalg.solve(S, a.mean(0) - b.mean(0))
    sa, sb = sc[pos], sc[~pos]
    return auroc(sc, pos), float(abs(sa.mean() - sb.mean()) / np.sqrt(0.5 * (sa.var() + sb.var())))


def pc1(X):
    Xc = X - X.mean(0); _, _, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Xc @ Vt[0]


def ridge_crosshalf(X, y, eo):
    """the supervised ceiling: y regressed on X (ridge) on one half, predicted on the other and vice versa; Spearman(y, prediction)."""
    pred = np.zeros(len(y))
    for fs, es in ((eo == 0, eo == 1), (eo == 1, eo == 0)):
        A = np.c_[X[fs], np.ones(fs.sum())]; lam = 1e-3 * np.trace(A.T @ A) / A.shape[1]
        w = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ y[fs]); pred[es] = np.c_[X[es], np.ones(es.sum())] @ w
    return spearman(y, pred)


def canon_corr(A, B):
    """canonical correlations between two coordinate sets (singular values of the whitened cross-covariance)."""
    A = A - A.mean(0); B = B - B.mean(0)
    Ca = np.linalg.cholesky(np.cov(A.T) + 1e-9 * np.eye(A.shape[1])); Cb = np.linalg.cholesky(np.cov(B.T) + 1e-9 * np.eye(B.shape[1]))
    M = np.linalg.solve(Ca, np.linalg.solve(Cb, (A.T @ B) / (len(A) - 1)).T).T
    return np.linalg.svd(M, compute_uv=False)


# ---------------------------------------------------------------- 1-D extreme deconvolution

def xd1(t, s2, K, maxit=300):
    """tied-variance K-mixture deconvolved from measurements t with variances s2 (latent_states_v4's XD on (t/s2, 1/s2))."""
    b = (t / s2)[:, None]; Gm = (1. / s2)[:, None, None]
    q = np.quantile(t, (np.arange(K) + 0.5) / K)
    r0 = LS.onehot(np.argmin(np.abs(t[:, None] - q[None]), 1), K) * 0.9 + 0.1 / K
    Om0 = np.array([[max(t.var() - s2.mean(), 0.05 * t.var())]]) / max(K, 1)
    return LS.xd_fit(b, Gm, r0, Om0, maxit=maxit)


def xd1_score(f, t, s2):
    return LS.xd_estep((t / s2)[:, None], (1. / s2)[:, None, None], f['pi'], f['xi'], f['Om'])[3]


def select_k(t, s2, eo, kmax):
    """cross-half held-out gain over K = 1 (both directions); K* = the K with the largest gain beyond 2 SE, else 1."""
    held = {}
    for K in range(1, kmax + 1):
        tot = np.zeros(0)
        for fs, es in ((eo == 0, eo == 1), (eo == 1, eo == 0)):
            f = xd1(t[fs], s2[fs], K); tot = np.r_[tot, xd1_score(f, t[es], s2[es])]
        held[K] = tot
    gains = {K: float((held[K] - held[1]).sum()) for K in held}
    ses = {K: float((held[K] - held[1]).std() * np.sqrt(len(held[K]))) for K in held}
    ks = [K for K in held if K > 1 and gains[K] - 2. * ses[K] > 0]
    return (max(ks, key=lambda K: gains[K] - 2. * ses[K]) if ks else 1), gains, ses


def mix_density(pi, xi, var, grid):
    return (pi[None] * np.exp(-0.5 * (grid[:, None] - xi[None]) ** 2 / var[None]) / np.sqrt(2 * np.pi * var[None])).sum(1)


def gap_depth(pi, xi, var, t, ga, gb):
    """D = 1 - valley / lower peak of the density between the children's modes: each child's mode is the highest local maximum of
    the deconvolved density inside its 5-95 % range (else the local maximum nearest its median); one shared mode means D = 0."""
    lo, hi = np.quantile(t, [0.005, 0.995]); grid = np.linspace(lo, hi, 600); fd = mix_density(pi, xi, var, grid)
    lm = np.where((fd[1:-1] > fd[:-2]) & (fd[1:-1] >= fd[2:]))[0] + 1
    if len(lm) == 0: lm = np.array([int(np.argmax(fd))])
    def mode(g):
        a, b = np.quantile(t[g], [0.05, 0.95]); inr = lm[(grid[lm] >= a) & (grid[lm] <= b)]
        j = inr[np.argmax(fd[inr])] if len(inr) else lm[np.argmin(np.abs(grid[lm] - np.median(t[g])))]
        return j
    ja, jb = mode(ga), mode(gb)
    if ja == jb: return 0., float(grid[ja]), float(grid[jb])
    if ja > jb: ja, jb = jb, ja
    valley = fd[ja:jb + 1].min(); pa, pb = fd[ja], fd[jb]
    return float(1. - valley / max(min(pa, pb), 1e-300)), float(grid[ja]), float(grid[jb])


def gap_from_fit(f, t, ga, gb):
    return gap_depth(f['pi'], f['xi'][:, 0], np.full(len(f['pi']), float(f['Om'][0, 0])), t, ga, gb)


def null_rep(args):
    """one uniform-continuum null replicate: latent ~ U(mu +- sqrt(3) sd), the observed noise variances permuted, same K and D."""
    seed, n, mu, sd, s2, K = args
    rng = np.random.default_rng(seed)
    lat = rng.uniform(mu - np.sqrt(3.) * sd, mu + np.sqrt(3.) * sd, n); s2p = rng.permutation(s2)
    t = lat + rng.normal(0., np.sqrt(s2p))
    return gap_from_fit(xd1(t, s2p, K, maxit=150), t, lat < mu, lat >= mu)[0]


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('tree'); ap.add_argument('out'); ap.add_argument('--truth', default='none')
    ap.add_argument('--nnull', type=int, default=100); ap.add_argument('--nproc', type=int, default=8); ap.add_argument('--kmax', type=int, default=8)
    ap.add_argument('--min_group', type=int, default=200); ap.add_argument('--max_pairs', type=int, default=12); ap.add_argument('--zsrc', default='ml')
    a = ap.parse_args(); T = a.tree.rstrip('/'); os.makedirs(a.out, exist_ok=True); t0 = time.time()
    rep = open(os.path.join(a.out, 'charts.txt'), 'w')
    def P(s): rep.write(s + '\n'); rep.flush(); log(s)
    nodes = [l.split()[0] for l in open(os.path.join(T, 'hetero_tree_nodes.txt')) if l.strip() and not l.startswith('#')]
    leaves = [l.split() for l in open(os.path.join(T, 'hetero_tree_leaves.txt')) if l.strip()]
    lab = read_ints(os.path.join(T, 'hetero_tree_labels.txt')); N = len(lab)
    leaf_node = {j + 1: w[0].split('_')[0] for j, w in enumerate(leaves)}
    kind, Tr = load_truth(a.truth, N)
    ch = {}
    P('# zsrc %s: per node the frame check (canonical correlations raw-MAP vs delivered z; |R^T R - I| of the LS rotation) and the variance budget per axis:' % a.zsrc)
    P('#   observed / flex deconvolved / measurement(ns P^-1 in the delivered frame) / [ml: un-shrunk observed - raw measurement]')
    for x in nodes:
        d = os.path.join(T, 'node_' + x); emb = os.path.join(d, 'fflex', 'flex_pca_embedding.bin'); co = os.path.join(d, 'fflex', 'flex_pca_coordinates.txt')
        if not (os.path.exists(emb) and os.path.exists(co)):
            P('node %s: no chart (fflex incomplete), skipped' % x); continue
        rows = read_ints(os.path.join(d, 'rows.txt')); E = read_embedding(emb); C = read_coords(co)
        assert E['n'] == len(rows) == len(C['z']) and (E['pinds'] == C['pind']).all(), 'node %s: rows / embedding / coordinates mismatch' % x
        k = E['k']; zt = C['z'][:, :k]; zb = E['z']
        pri = read_deconv_prior(os.path.join(d, 'fflex', 'flex_pca_deconv_prior.txt')) if os.path.exists(os.path.join(d, 'fflex', 'flex_pca_deconv_prior.txt')) else None
        ns = pri['noise_scale'] if pri is not None else E['noise_scale']
        # the delivered frame = the raw one rotated (gauge fix): least-squares rotation and its orthogonality
        X = np.c_[zb, np.ones(len(zb))]; Bm = np.linalg.lstsq(X, zt, rcond=None)[0]; R = Bm[:k]; orth = float(np.abs(R.T @ R - np.eye(k)).max())
        cc = canon_corr(zt, zb); resid = float(np.sqrt(((zt - X @ Bm) ** 2).sum(1).mean()) / np.sqrt((zt - zt.mean(0)).var(0).sum()))
        cov_raw = np.linalg.inv(E['prec'])                                      # likelihood covariance of the raw MAP embedding
        # the budget that balances (observed = flex deconvolved + measurement) is the un-shrunk z with the noise-scaled covariance
        if a.zsrc == 'ml':
            pp = read_map_prior(os.path.join(d, 'fflex', 'flex_pca_map_prior.txt'))      # prior precision per component
            zml = np.einsum('nij,nj->ni', cov_raw @ (E['prec'] + np.diag(pp[:k])[None]), zb)    # z_ML = P^-1 (P + Pi) z_MAP
            z = np.c_[zml, np.ones(len(zml))] @ Bm
        else:
            z = zt
        cov = np.einsum('ji,njk,kl->nil', R, cov_raw, R) * ns
        vo = z.var(0); vm = np.array([cov[:, j, j].mean() for j in range(k)])
        vd = (pri['w'][:, None] * (pri['sig'] + pri['mu'] ** 2)).sum(0) - ((pri['w'][:, None] * pri['mu']).sum(0)) ** 2 if pri is not None else np.full(k, np.nan)
        pre = read_ints(os.path.join(d, 'split', 'split_children_pre.txt')); acc = read_ints(os.path.join(d, 'split', 'split_children.txt'))
        ch[x] = dict(rows=rows, z=z, cov=cov, eo=C['eo'], pri=pri, pre=pre, acc=acc, loc={r: i for i, r in enumerate(rows)}, k=k, ns=ns, R=R, Bm=Bm)
        P('node %-6s n %6d k %2d ns %.3f | canon corr %s, rotation residual %.3f, |R^T R - I| %.2f | axes 1-5: %s' % (
            x, len(rows), k, ns, ' '.join('%.3f' % v for v in cc[:4]), resid, orth, ' '.join('%.0f/%.0f/%.0f' % (vo[j], vd[j], vm[j]) for j in range(min(k, 5)))))
    summary = dict(tree=T, kind=kind, zsrc=a.zsrc, nodes=list(ch), coordination={}, splits=[])

    # ---------------- 1. coordination
    maps = {}
    for x in ch:
        if x == 'r' or x[:-1] not in ch: continue
        p = x[:-1]; c = ch[x]; pc = ch[p]
        ip = np.array([pc['loc'][r] for r in c['rows']]); Zp = pc['z'][ip]; Zc = c['z']
        w = 1. / np.array([np.trace(pc['cov'][i]) / pc['k'] for i in ip]); w /= w.mean()
        X = np.c_[Zc, np.ones(len(Zc))]
        for tag, ww in (('w', w), ('u', np.ones(len(Zc)))):
            B = np.linalg.solve((X * ww[:, None]).T @ X + 1e-9 * np.eye(X.shape[1]), (X * ww[:, None]).T @ Zp)
            maps[(x, tag)] = B
            if tag == 'w':
                r2 = 1. - (ww[:, None] * (Zp - X @ B) ** 2).sum(0) / (ww[:, None] * (Zp - Zp.mean(0)) ** 2).sum(0)
                P('map %s -> %s: %d particles, R^2 per parent axis %s' % (x, p, len(Zc), ' '.join('%.2f' % v for v in r2[:6])))
    def compose(x, Z, tag):
        while x != 'r':
            Z = np.c_[Z, np.ones(len(Z))] @ maps[(x, tag)]; x = x[:-1]
        return Z
    deepest = np.array([leaf_node.get(l, 'r') if l > 0 else 'r' for l in lab], dtype=object)
    for i in range(N):
        if deepest[i] not in ch: deepest[i] = 'r'
    rootloc = ch['r']['loc']; have = np.array([r in rootloc for r in range(1, N + 1)])
    G = {}
    for tag in ('w', 'u'):
        Z = np.full((N, ch['r']['k']), np.nan)
        for x in ch:
            sel = np.where((deepest == x) & have)[0]
            if len(sel) == 0: continue
            loc = np.array([ch[x]['loc'][r + 1] for r in sel]); Z[sel] = compose(x, ch[x]['z'][loc], tag)
        G['global_' + tag] = Z
    Zr = np.full((N, ch['r']['k']), np.nan); ir = np.where(have)[0]; Zr[ir] = ch['r']['z'][[rootloc[r + 1] for r in ir]]; G['root'] = Zr
    eo_all = np.full(N, -1); eo_all[ir] = ch['r']['eo'][[rootloc[r + 1] for r in ir]]
    if kind == 'qtoy':
        arm = Tr['arm'] & have; closed = (Tr['conformer'] == 18) & have; holo = (Tr['holo'] > 0.5) & have
        conf = Tr['conformer'].astype(int); s1 = Tr['s1']
        P('\n== coordination (%d arm, %d closed, %d holo particles): |Spearman(s1, PC1 of the arm particles)| overall / per conformer mean, the supervised s1 ceiling (cross-half ridge, Spearman); closed-vs-arm and holo-vs-arm cross-half LDA AUROC (d\')' % (arm.sum(), closed.sum(), holo.sum()))
        for name, Z in G.items():
            ok = np.isfinite(Z).all(1)
            pa = pc1(Z[arm & ok]); rho = spearman(s1[arm & ok], pa); ceil = ridge_crosshalf(Z[arm & ok], s1[arm & ok], eo_all[arm & ok])
            cm = {c: pa[conf[arm & ok] == c].mean() for c in np.unique(conf[arm & ok])}
            rho_c = spearman(np.array([s1[conf == c][0] for c in cm]), np.array(list(cm.values()))) if len(cm) > 2 else float('nan')
            ca, cd = lda_crosshalf(Z[(arm | closed) & ok], closed[(arm | closed) & ok], eo_all[(arm | closed) & ok])
            ha, hd = lda_crosshalf(Z[(arm | holo) & ok], holo[(arm | holo) & ok], eo_all[(arm | holo) & ok])
            summary['coordination'][name] = dict(spearman_s1=rho, spearman_conformer_means=rho_c, s1_ceiling=ceil, closed_auroc=ca, closed_dprime=cd, holo_auroc=ha, holo_dprime=hd)
            P('  %-9s rho(s1) %.3f  rho(conformer means) %.3f  ceiling %.3f | closed AUROC %.3f d\' %.2f | holo AUROC %.3f d\' %.2f' % (name, abs(rho), abs(rho_c), ceil, ca, cd, ha, hd))
        P('  per node (arm particles of the node): |rho(s1)| of PC1 local chart / composed(w) / root-direct, then the s1 ceiling local / composed / root; n')
        for x in ch:
            sel = np.where((deepest == x) & arm)[0]
            if len(sel) < 50: continue
            loc = np.array([ch[x]['loc'][r + 1] for r in sel]); zl = ch[x]['z'][loc]; eo_s = eo_all[sel]
            rl = spearman(s1[sel], pc1(zl)); rc = spearman(s1[sel], pc1(G['global_w'][sel])); rr = spearman(s1[sel], pc1(G['root'][sel]))
            cl = ridge_crosshalf(zl, s1[sel], eo_s); cc = ridge_crosshalf(G['global_w'][sel], s1[sel], eo_s); cr = ridge_crosshalf(G['root'][sel], s1[sel], eo_s)
            P('    %-6s PC1 %.3f / %.3f / %.3f   ceiling %.3f / %.3f / %.3f  (%d; s1 5-95 %% %.1f .. %.1f)' % (x, abs(rl), abs(rc), abs(rr), cl, cc, cr, len(sel), np.percentile(s1[sel], 5), np.percentile(s1[sel], 95)))
            summary['coordination']['node_' + x] = dict(local=rl, composed=rc, root=rr, ceiling_local=cl, ceiling_composed=cc, ceiling_root=cr, n=int(len(sel)))

    # ---------------- 2. split criterion
    P('\n== split criterion per node and child pair: separation Delta/sd(obs) and d\' (deconvolved, floored at 10 %% of the observed variance), signal fraction of the split direction, K*, D (pooled; even / odd), p vs the uniform-continuum null (%d replicates), D_flex from flex\'s deconvolved prior; composition' % a.nnull)
    pool = Pool(a.nproc)
    for x in ch:
        c = ch[x]; g = c['pre'] if c['pre'] is not None else c['acc']
        if g is None: continue
        groups = [u for u in np.unique(g) if u > 0 and (g == u).sum() >= a.min_group]
        if len(groups) < 2: continue
        accepted = c['acc'] is not None and len(np.unique(c['acc'][c['acc'] > 0])) > 1
        for ga, gb in [(ga, gb) for i, ga in enumerate(groups) for gb in groups[i + 1:]][:a.max_pairs]:
            ma, mb = g == ga, g == gb; sel = ma | mb
            u = c['z'][mb].mean(0) - c['z'][ma].mean(0); u /= max(np.linalg.norm(u), 1e-12)
            t = c['z'] @ u; s2 = np.maximum(np.einsum('i,nij,j->n', u, c['cov'], u), 1e-6 * t.var())
            ts, s2s, eos = t[sel], s2[sel], c['eo'][sel]
            delta = t[mb].mean() - t[ma].mean(); sdo = np.sqrt(0.5 * (t[ma].var() + t[mb].var()))
            vdec = 0.5 * sum(max(t[m].var() - s2[m].mean(), 0.1 * t[m].var()) for m in (ma, mb)); dprime = abs(delta) / np.sqrt(vdec)
            sfrac = float(max(0., 1. - s2s.mean() / ts.var()))
            kstar, gains, ses = select_k(ts, s2s, eos, a.kmax)
            f = xd1(ts, s2s, kstar); D, moda, modb = gap_from_fit(f, ts, ma[sel], mb[sel])
            Dh = [gap_from_fit(xd1(ts[eos == h], s2s[eos == h], kstar), ts[eos == h], ma[sel][eos == h], mb[sel][eos == h])[0] for h in (0, 1)]
            mu = float((f['pi'] * f['xi'][:, 0]).sum()); sd = float(np.sqrt(max((f['pi'] * (f['xi'][:, 0] ** 2 + f['Om'][0, 0])).sum() - mu ** 2, 1e-12)))
            Dn = np.array(pool.map(null_rep, [(1000 + j, len(ts), mu, sd, s2s, kstar) for j in range(a.nnull)]))
            pval = (1. + (Dn >= D).sum()) / (a.nnull + 1.)
            Df = float('nan')
            if c['pri'] is not None:                                             # flex's own deconvolved prior on the same direction (delivered frame)
                uu = u if a.zsrc != 'ml' else u
                Df = gap_depth(c['pri']['w'], c['pri']['mu'][:, :c['k']] @ uu, np.maximum((c['pri']['sig'][:, :c['k']] * uu ** 2).sum(1), 1e-12), ts, ma[sel], mb[sel])[0]
            comp = ''
            if kind == 'qtoy':
                def desc(m):
                    rr = c['rows'][m] - 1; cf = Tr['conformer'].astype(int)[rr]; jk = Tr['junk'][rr] > 0.5
                    u_, n_ = np.unique(cf[~jk], return_counts=True); o = np.argsort(-n_)[:2]
                    s1m = Tr['s1'][rr][Tr['arm'][rr]]
                    return ' '.join('%d(%.2f)' % (u_[i], n_[i] / max((~jk).sum(), 1)) for i in o) + ' junk %.2f' % jk.mean() + (' s1 %.1f..%.1f' % (np.percentile(s1m, 10), np.percentile(s1m, 90)) if len(s1m) > 20 else '')
                comp = ' | a: %s; b: %s' % (desc(ma), desc(mb))
            elif kind == 'closed':
                ref = Tr['closedref'][c['rows'] - 1] > 0.5
                comp = ' | closed refs a %d/%d b %d/%d (enrichment %.2f / %.2f)' % (ref[ma].sum(), ma.sum(), ref[mb].sum(), mb.sum(), ref[ma].mean() / max(ref.mean(), 1e-9), ref[mb].mean() / max(ref.mean(), 1e-9))
            P('  %-5s pair %d|%d (%s) n %d %d: Delta/sd %.2f d\' %.2f sfrac %.2f K* %d D %.3f (%.3f %.3f) null med %.3f p95 %.3f p %.3f D_flex %.3f%s' % (
                x, ga, gb, 'accepted' if accepted else 'merged', ma.sum(), mb.sum(), abs(delta) / sdo, dprime, sfrac, kstar, D, Dh[0], Dh[1], np.median(Dn), np.percentile(Dn, 95), pval, Df, comp))
            summary['splits'].append(dict(node=x, pair=[int(ga), int(gb)], accepted=bool(accepted), n=[int(ma.sum()), int(mb.sum())], delta_sd=float(abs(delta) / sdo), dprime=float(dprime),
                                          signal_fraction=sfrac, kstar=int(kstar), gains={int(k_): v for k_, v in gains.items()}, D=D, D_halves=Dh, null_median=float(np.median(Dn)),
                                          null_p95=float(np.percentile(Dn, 95)), p=float(pval), D_flex=Df, composition=comp))
    pool.close(); pool.join()
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'charts_summary.json'), 'w'), indent=1, default=float)
    P('done in %.0f s' % summary['wall']); rep.close()


if __name__ == '__main__':
    main()

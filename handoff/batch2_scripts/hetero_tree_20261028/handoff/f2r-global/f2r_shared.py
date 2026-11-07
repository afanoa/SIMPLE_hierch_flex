#!/usr/bin/env python3
"""f2r_shared.py: a seam-free global coordinate. ONE tree-wide basis: every charted node's flex modes (band-limited to flex's 15 A),
orthonormalised in the band, the top q directions by SVD of their Gram matrix (--basis root keeps the root's modes only: the
root-chart control in the same pipeline). Every particle is embedded at its own pose (the root project's input poses, or the
final leaf stage-B poses assembled by f2r_leafposes.py) with its own sufficient statistics b, G (latent_states_v4): z_ML with
covariance (G + ridge)^-1, and z_MAP under one Gaussian fitted to all particles. Tests as f2r_global.py: (a) closed AUROC / d',
per-leaf isolation, a named pair; (c) leaf pairs vs truth equivalence; (d) valleys on the root split axis and on the named pair's
axis; (e) the Bhattacharyya merge criterion over the reported pairs.
usage: f2r_shared.py embed --tree T --node T/node_r --dump D --cache C --out O [--truth K:ARGS --basis all|root --q 30 --poses input|final
       --pairs all|siblings --merged F --pair A,B --dsub N --maxleaves 12]
       f2r_shared.py compare OUT_A OUT_B NAME_A NAME_B"""
import argparse, glob, json, os, sys, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.environ.get('MEST_DIR', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest'))
import latent_states_v4 as LS
import flexprior_estep as FP
import f2r_calib as FC
import f2r_global as FG
import mest_charts as MC

G = LS.G
log = LS.log


def tree_basis(T, nodes, kf, which, qmax, P):
    """the band-orthonormal top-q span of the nodes' flex modes (all charted nodes, or the root alone)."""
    from mrcio import read_mrc
    band = np.r_[(G['cell_shell'] <= G['kb']), False]
    fw = np.r_[G['cell_fw'], 0.] * band
    rows, src = [], []
    for x in nodes:
        if which == 'root' and x != 'r': continue
        for f in sorted(glob.glob(os.path.join(T, 'node_' + x, 'fflex', 'flex_pca_pc[0-9][0-9][0-9].mrc'))):
            V = FP.lowpass_cells(LS.vol_ft(read_mrc(f)[0]), kf)
            rows.append(V / np.sqrt(max(np.sum(fw * np.abs(V) ** 2), 1e-300))); src.append(x)
    A = np.array(rows); nn = len(rows)
    As = A * np.sqrt(fw)[None]
    M = (As @ As.conj().T).real
    del As
    lam, U = np.linalg.eigh(M); order = np.argsort(-lam); lam, U = np.maximum(lam[order], 0.), U[:, order]
    keep = int((lam > 1e-6 * lam[0]).sum()); q = min(qmax, keep)
    W = (U[:, :q].T @ A) / np.sqrt(lam[:q])[:, None]
    cum = np.cumsum(lam) / lam.sum()
    qat = [int(np.searchsorted(cum, c) + 1) for c in (0.9, 0.95, 0.99)]
    # how much of each root mode the kept span holds
    ir = [i for i, s in enumerate(src) if s == 'r']
    cover = [float(np.sum(np.abs((W.conj() * fw[None]) @ A[i]) ** 2)) for i in ir]
    P('basis %s: %d maps from %d nodes, rank %d -> q %d holds %.3f of the energy (q at 90/95/99 %%: %d/%d/%d); the root modes'' share inside the span: %s'
      % (which, nn, len(set(src)), keep, q, cum[q - 1], qat[0], qat[1], qat[2], ' '.join('%.2f' % c for c in cover)))
    return [W[l] for l in range(q)], dict(nmaps=nn, nnodes=len(set(src)), rank=keep, q=q, energy=float(cum[q - 1]), q90=qat[0], q95=qat[1], q99=qat[2], root_cover=cover)


def embed_cmd(a):
    T = a.tree.rstrip('/'); os.makedirs(a.out, exist_ok=True); t0 = time.time()
    rep = open(os.path.join(a.out, 'report.txt'), 'w')
    def P(s): rep.write(s + '\n'); rep.flush(); log(s)
    a.lp = FP.LPB; a.poses = 'stageB' if a.poses == 'final' else 'input'; a.mskdiam = 0.; a.maxn = 0
    FC._POOL[:] = [ThreadPoolExecutor(a.nproc), a.nproc]
    LS.xd_estep, LS.xd_mstep = FC.xd_estep_par, FC.xd_mstep_fast
    LS.setup(a)
    from mrcio import read_mrc
    from ppinfo_site import soft_sphere
    n, n0, didx, prow = G['n'], G['n0'], G['didx'], G['prow']
    G['sm_c'] = np.fft.ifftshift(soft_sphere(G['B'], G['msk_crop']))
    kf = int(round(G['box_A'] / FP.LPF))
    nodes = [l.split()[0] for l in open(os.path.join(T, 'hetero_tree_nodes.txt')) if l.strip() and not l.startswith('#')]
    charted = [x for x in nodes if os.path.exists(os.path.join(T, 'node_' + x, 'fflex', 'flex_pca_pc001.mrc'))]
    leaves = [l.split()[0] for l in open(os.path.join(T, 'hetero_tree_leaves.txt')) if l.strip()]
    lab_all = MC.read_ints(os.path.join(T, 'hetero_tree_labels.txt')); assert len(lab_all) == n0
    lab = lab_all[prow - 1]; names = leaves; L = len(names)
    P('tree %s: %d particles embedded (project %d rows), %d leaves, %d charted nodes, poses %s, flex band to ring %d of %d'
      % (T, n, n0, L, len(charted), a.poses, kf, G['kb']))
    W, binfo = tree_basis(T, charted, kf, a.basis, a.q, P)
    V0 = LS.vol_ft(read_mrc(G['vol1'])[0])
    beta_mu, sig2_mu, amp = LS.scale_noise(V0)
    G['cn'] = amp
    b, Gm = LS.suff_stats(V0, W, beta_mu, sig2_mu, amp)
    G.update(b=b, Gm=Gm)
    q = b.shape[1]
    ridge = 1e-3 * np.median(np.trace(Gm, axis1=1, axis2=2)) / q
    Gr = Gm + ridge * np.eye(q)[None]
    cov = np.linalg.inv(Gr); z = np.einsum('nij,nj->ni', cov, b)
    G['mu1'] = np.zeros(q)
    f1 = LS.fit_one(np.arange(n), np.eye(q) / max(np.median(np.trace(Gm, axis1=1, axis2=2)) / q, 1e-30) * 10.)
    Om1, mu1 = f1['Om'], f1['xi'][0]; Oi = np.linalg.inv(Om1)
    zmap = np.einsum('nij,nj->ni', np.linalg.inv(Gm + Oi[None]), b + Oi @ mu1)
    vo = z.var(0); vm = np.array([cov[:, j, j].mean() for j in range(q)])
    P('embedding: q %d; observed / measurement variance, axes 1-6: %s; one-Gaussian LL/particle %.3f' % (q, ' '.join('%.0f/%.0f' % (vo[j], vm[j]) for j in range(min(q, 6))), f1['LL'] / n))
    eo = G['eo']
    kind, Tr = MC.load_truth(a.truth, n0)
    ok = np.ones(n, bool)
    truth_cols = {}
    closed = None
    if kind == 'closed':
        closed = Tr['closedref'][prow - 1] > 0.5; truth_cols['closed'] = closed
    elif kind == 'qtoy':
        truth_cols = {k_: Tr[k_][prow - 1] for k_ in ('conformer', 'combo', 'arm', 's1')}
        closed = truth_cols['conformer'].astype(int) == 18
        if closed.sum() == 0: closed = None
    summary = dict(tree=T, basis=a.basis, poses=a.poses, n=int(n), q=int(q), kind=kind, basis_info=binfo)
    # the pairs to report
    allowed = None; merged_pairs = set(); named = None
    if a.pairs == 'siblings' or a.merged or a.pair:
        allowed = {(i, j) for i in range(1, L + 1) for j in range(i + 1, L + 1) if names[i - 1].split('_')[0] == names[j - 1].split('_')[0]}
        if a.pairs == 'all': allowed = {(i, j) for i in range(1, L + 1) for j in range(i + 1, L + 1)}
    if a.merged:
        mu = MC.read_ints(a.merged)[prow - 1]; unit = {}
        for i in range(1, L + 1):
            m = (lab == i) & (mu > 0)
            if m.sum(): u, c = np.unique(mu[m], return_counts=True); unit[i] = int(u[np.argmax(c)])
        merged_pairs = {(i, j) for i in unit for j in unit if i < j and unit[i] == unit[j]}; allowed |= merged_pairs
    if a.pair:
        n1, n2 = a.pair.split(','); i, j = sorted((names.index(n1) + 1, names.index(n2) + 1)); named = (i, j); allowed.add(named)
    for zname, Z in (('z_ML', z), ('z_MAP', zmap)):
        if closed is not None:
            au, dp = MC.lda_crosshalf(Z, closed, eo)
            summary['closed_' + zname] = dict(auroc=au, dprime=dp)
            P('== (a) %s: closed (%d of %d) vs the rest, cross-half LDA AUROC %.3f d\' %.2f' % (zname, closed.sum(), n, au, dp))
    # per-leaf isolation and pairs on z_ML
    pg = FG.leaf_pairs(z, eo, lab, names, lambda s: None, 'shared', allowed=allowed)
    if closed is not None:
        share = closed.mean(); rows = []
        order = sorted(range(1, L + 1), key=lambda i: -closed[lab == i].mean() if (lab == i).sum() else 0.)
        P('   per leaf: n, closed count, enrichment; min d\' to another reported leaf')
        for i in order:
            m = lab == i
            if m.sum() == 0: continue
            dg = [pg[k_]['dprime'] for k_ in pg if i in k_]
            rows.append(dict(leaf=names[i - 1], n=int(m.sum()), nclosed=int(closed[m].sum()), enrichment=float(closed[m].mean() / max(share, 1e-12)), min_dprime=float(min(dg)) if dg else None))
            if len(rows) <= a.maxleaves:
                P('   %-12s n %7d closed %5d enrichment %5.2f | min d\' %s' % (names[i - 1], m.sum(), closed[m].sum(), closed[m].mean() / max(share, 1e-12), '%.2f' % min(dg) if dg else '-'))
        summary['leaves'] = rows
    if L <= a.maxleaves:
        P('   leaf pairs (shared frame, z_ML):')
        FG.leaf_pairs(z, eo, lab, names, P, 'shared', allowed=allowed)
    ds = sorted(v['dprime'] for v in pg.values())
    if ds:
        P('   %d reported pairs: d\' quartiles %.2f / %.2f / %.2f, min %.2f max %.2f' % (len(ds), np.percentile(ds, 25), np.median(ds), np.percentile(ds, 75), ds[0], ds[-1]))
    if merged_pairs:
        dm = sorted(pg[k_]['dprime'] for k_ in pg if k_ in merged_pairs); dn = sorted(pg[k_]['dprime'] for k_ in pg if k_ not in merged_pairs)
        P('   pairs united by the merge file: d\' %s' % ' '.join('%.2f' % v for v in dm))
        if dn: P('   other reported pairs: d\' quartiles %.2f / %.2f / %.2f (%d)' % (np.percentile(dn, 25), np.median(dn), np.percentile(dn, 75), len(dn)))
        summary['merged_pairs_dprime'] = dm; summary['other_pairs_dprime'] = dn
    if named is not None and named in pg:
        v = pg[named]; rank = sum(1 for x in pg.values() if x['dprime'] > v['dprime']); rbc = sum(1 for x in pg.values() if x['bc'] < v['bc'])
        P('   named pair %s: d\' %.2f AUROC %.3f overlap %.3f; %d of %d reported pairs further apart by d\', %d less overlapping; united by the merge file: %s'
          % (a.pair, v['dprime'], v['auroc'], v['bc'], rank, len(pg), rbc, named in merged_pairs))
        summary['named_pair'] = dict(pair=a.pair, **v, rank=rank, rank_bc=rbc, npairs=len(pg), merged=named in merged_pairs)
    # (c) truth equivalence
    if kind == 'qtoy':
        key = 'combo' if len(np.unique(truth_cols['combo'])) <= 6 else 'conformer'; tv = truth_cols[key].astype(int); dom = {}
        for i in range(1, L + 1):
            m = lab == i
            if m.sum() == 0: continue
            u, c = np.unique(tv[m], return_counts=True); dom[i] = (int(u[np.argmax(c)]), float(c.max() / m.sum()))
        equiv = {(i, j): dom[i][0] == dom[j][0] for i in dom for j in dom if i < j}
        de = [v['dprime'] for k_, v in pg.items() if equiv.get(k_)]; dn = [v['dprime'] for k_, v in pg.items() if equiv.get(k_) is False]
        P('== (c) dominant %s per leaf: %s' % (key, ' '.join('%s:%d(%.2f)' % (names[i - 1], d[0], d[1]) for i, d in dom.items())))
        if de and dn:
            sc = np.array(de + dn); pos = np.r_[np.zeros(len(de), bool), np.ones(len(dn), bool)]
            P('   d\' of equivalent pairs %s | distinct pairs min / median %.2f / %.2f | AUROC(d\' separates distinct from equivalent) %.3f'
              % (' '.join('%.2f' % v for v in sorted(de)), min(dn), np.median(dn), MC.auroc(sc, pos)))
            summary['equiv'] = dict(equivalent=de, distinct=dn, auroc=MC.auroc(sc, pos))
        else:
            P('   no truth-equivalent leaf pair in this tree' if not de else '   no distinct pair')
    # (d) valleys: the root split axis, and the named pair's axis
    acc = MC.read_ints(os.path.join(T, 'node_r', 'split', 'split_children.txt'))
    rs = np.zeros(n, int)
    if acc is not None and len(acc) == n: rs = acc.copy()
    elif acc is not None and len(acc) == n0: rs = acc[prow - 1]
    groups = [g for g in np.unique(rs) if g > 0]
    subset = truth_cols['arm'] if 'arm' in truth_cols else np.ones(n, bool)
    if a.dsub and subset.sum() > a.dsub:
        keep = np.zeros(n, bool); keep[np.random.default_rng(0).choice(np.where(subset)[0], a.dsub, replace=False)] = True; subset = subset & keep
    def axis_test(ma, mb, sel, label, cutlabel):
        d = z[ma].mean(0) - z[mb].mean(0); d /= np.linalg.norm(d)
        t = z[sel] @ d; s2 = np.einsum('i,nij,j->n', d, cov[sel], d); s2 = np.maximum(s2, 1e-6 * t.var())
        cut = 0.5 * (np.median(z[ma] @ d) + np.median(z[mb] @ d))
        P('== (d) %s: %d / %d particles, the %s at %+.2f sd; deconvolved density per half:' % (label, ma.sum(), mb.sum(), cutlabel, (cut - t.mean()) / t.std()))
        return dict(cut_sd=float((cut - t.mean()) / t.std()), valleys=FG.valleys(t, s2, eo[sel], a.kmax, P, label))
    if len(groups) >= 2:
        sizes = {g: int((rs == g).sum()) for g in groups}; g1, g2 = sorted(groups, key=lambda g: -sizes[g])[:2]
        summary['cut_root'] = axis_test(subset & (rs == g1), subset & (rs == g2), subset, 'root split %d|%d axis' % (g1, g2), "tree's cut")
    if named is not None:
        i, j = named; summary['cut_pair'] = axis_test(lab == i, lab == j, (lab == i) | (lab == j), 'named pair axis', 'midpoint')
    # (e) merge criterion
    if pg:
        bcs = sorted(v['bc'] for v in pg.values())
        P('== (e) Bhattacharyya overlap of the reported leaf pairs: quartiles %.3f / %.3f / %.3f; pairs with overlap > 0.5: %d of %d%s'
          % (np.percentile(bcs, 25), np.median(bcs), np.percentile(bcs, 75), sum(v > 0.5 for v in bcs), len(bcs),
             ('; united-by-merge pairs > 0.5: %d of %d' % (sum(pg[k_]['bc'] > 0.5 for k_ in merged_pairs if k_ in pg), len([k_ for k_ in merged_pairs if k_ in pg]))) if merged_pairs else ''))
    np.savez(os.path.join(a.out, 'shared.npz'), z=z, zmap=zmap, covdiag=np.einsum('nii->ni', cov), lab=lab, prow=prow, eo=eo, names=np.array(names), rootsplit=rs,
             **{k_: v for k_, v in truth_cols.items()})
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1, default=float)
    P('done in %.0f s' % summary['wall']); rep.close()
    FC._POOL[0].shutdown()


def compare_cmd(a):
    A = np.load(os.path.join(a.out_a, 'shared.npz'), allow_pickle=True); B = np.load(os.path.join(a.out_b, 'shared.npz'), allow_pickle=True)
    common, ia, ib = np.intersect1d(A['prow'], B['prow'], return_indices=True)
    la, lb = A['lab'][ia], B['lab'][ib]
    print('== %s vs %s: %d shared particles; leaf-partition NMI %.3f' % (a.name_a, a.name_b, len(common), LS.nmi(la, lb)))
    def whiten(X):
        w, V = np.linalg.eigh(np.cov(X.T)); keep = w > 1e-8 * w.max()
        return X @ V[:, keep] / np.sqrt(w[keep])
    for key in ('z', 'zmap'):
        Xa = A[key][ia]; Xb = B[key][ib]; Xa = Xa - Xa.mean(0); Xb = Xb - Xb.mean(0)
        Wa, Wb = whiten(Xa), whiten(Xb)
        U, S, Vt = np.linalg.svd(Wa.T @ Wb / (len(common) - 1), full_matrices=False)
        rhos = []; top = bot = float('nan')
        for c in range(min(3, len(S))):
            ua, ub = Wa @ U[:, c], Wb @ Vt[c]
            if np.corrcoef(ua, ub)[0, 1] < 0: ub = -ub
            rhos.append(MC.spearman(ua, ub))
            if c == 0:
                n5 = int(0.05 * len(common))
                top = len(set(np.argsort(ua)[-n5:]) & set(np.argsort(ub)[-n5:])) / n5; bot = len(set(np.argsort(ua)[:n5]) & set(np.argsort(ub)[:n5])) / n5
        print('   %-5s canonical correlations %s | variates 1-3 Spearman %s; variate 1 top-5%% overlap %.2f, bottom-5%% %.2f (chance 0.05)'
              % (key, ' '.join('%.3f' % v for v in S[:6]), ' '.join('%.3f' % r for r in rhos), top, bot))
    if 'closed' in A.files:
        cl = A['closed'][ia]
        ea = {i: cl[la == i].mean() / cl.mean() for i in np.unique(la) if i > 0}; eb = {i: cl[lb == i].mean() / cl.mean() for i in np.unique(lb) if i > 0}
        i_, j_ = max(ea, key=ea.get), max(eb, key=eb.get)
        inter = ((la == i_) & (lb == j_)).sum()
        print('   most closed-enriched leaves %s (x%.2f) / %s (x%.2f): %d shared particles of %d / %d (Jaccard %.2f)'
              % (A['names'][i_ - 1], ea[i_], B['names'][j_ - 1], eb[j_], inter, (la == i_).sum(), (lb == j_).sum(), inter / ((la == i_) | (lb == j_)).sum()))


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    e = sub.add_parser('embed')
    for k_ in ('tree', 'node', 'dump', 'out'): e.add_argument('--' + k_, required=True)
    e.add_argument('--cache', default=''); e.add_argument('--truth', default='none'); e.add_argument('--basis', default='all')
    e.add_argument('--q', type=int, default=30); e.add_argument('--poses', default='final'); e.add_argument('--nproc', type=int, default=8)
    e.add_argument('--pairs', default='all'); e.add_argument('--merged', default=''); e.add_argument('--pair', default='')
    e.add_argument('--dsub', type=int, default=0); e.add_argument('--maxleaves', type=int, default=12); e.add_argument('--kmax', type=int, default=6)
    c = sub.add_parser('compare'); c.add_argument('out_a'); c.add_argument('out_b'); c.add_argument('name_a'); c.add_argument('name_b')
    a = ap.parse_args()
    embed_cmd(a) if a.cmd == 'embed' else compare_cmd(a)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""f2r_global.py: tests of the global chart coordinate (mest_charts.py's coordination: every node's fflex fit is a local PPCA chart,
child -> parent precision-weighted affine maps composed to the root; a particle's global coordinate is its deepest chart's un-shrunk
z carried into the root frame).
  compose TREE OUT --truth KIND:ARGS   per-particle global and root-chart coordinates (OUT/global.npz) and, in OUT/report.txt:
     (a) minority: closed-vs-rest cross-half LDA AUROC / d' in the global frame and the root chart; per leaf the closed enrichment,
         and the leaf's cross-half d' to every other leaf (does the closed leaf stand apart?);
     (c) leaf pairs: cross-half d', AUROC and the Bhattacharyya overlap of Gaussian fits in the global frame, against the truth
         equivalence of the pair (same dominant conformer / kind), and the merge decisions a d' threshold would make;
     (d) the cut of the root split along the global axis: the deconvolved 1-D density (XD, K by cross-half held-out) of the arm
         particles on the split direction, fitted on each particle half, its valleys (depth, position) against the tree's cut.
  compare OUT_A OUT_B NAME_A NAME_B      (b) reproducibility: canonical correlations of the two trees' global (and root) coordinates
     on the shared particles, per-particle agreement of the first canonical variates, extrema (top / bottom 5 %) overlap, closed-leaf
     correspondence, against the leaf-partition NMI."""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.environ.get('MEST_DIR', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/mest'))
sys.path.insert(0, os.environ.get('MEST_LATENT', '/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent'))
import mest_charts as MC
import latent_states_v4 as LS

log = LS.log


# ---------------------------------------------------------------- the tree's charts and their composition (mest_charts.py, zsrc=ml)

def load_tree(T):
    nodes = [l.split()[0] for l in open(os.path.join(T, 'hetero_tree_nodes.txt')) if l.strip() and not l.startswith('#')]
    leaves = [l.split() for l in open(os.path.join(T, 'hetero_tree_leaves.txt')) if l.strip()]
    lab = MC.read_ints(os.path.join(T, 'hetero_tree_labels.txt')); N = len(lab)
    leaf_names = [w[0] for w in leaves]
    leaf_node = {j + 1: w[0].split('_')[0] for j, w in enumerate(leaves)}
    ch = {}
    for x in nodes:
        d = os.path.join(T, 'node_' + x); emb = os.path.join(d, 'fflex', 'flex_pca_embedding.bin'); co = os.path.join(d, 'fflex', 'flex_pca_coordinates.txt')
        if not (os.path.exists(emb) and os.path.exists(co)):
            log('node %s: no chart, skipped' % x); continue
        rows = MC.read_ints(os.path.join(d, 'rows.txt')); E = MC.read_embedding(emb); C = MC.read_coords(co)
        assert E['n'] == len(rows) == len(C['z']) and (E['pinds'] == C['pind']).all()
        k = E['k']; zt = C['z'][:, :k]; zb = E['z']
        prf = os.path.join(d, 'fflex', 'flex_pca_deconv_prior.txt')
        ns = MC.read_deconv_prior(prf)['noise_scale'] if os.path.exists(prf) else E['noise_scale']
        X = np.c_[zb, np.ones(len(zb))]; Bm = np.linalg.lstsq(X, zt, rcond=None)[0]; R = Bm[:k]
        cov_raw = np.linalg.inv(E['prec'])
        pp = MC.read_map_prior(os.path.join(d, 'fflex', 'flex_pca_map_prior.txt'))
        zml = np.einsum('nij,nj->ni', cov_raw @ (E['prec'] + np.diag(pp[:k])[None]), zb)
        z = np.c_[zml, np.ones(len(zml))] @ Bm
        cov = np.einsum('ji,njk,kl->nil', R, cov_raw, R) * ns
        ch[x] = dict(rows=rows, z=z, cov=cov, eo=C['eo'], loc={r: i for i, r in enumerate(rows)}, k=k, ns=ns,
                     acc=MC.read_ints(os.path.join(d, 'split', 'split_children.txt')))
    maps = {}
    for x in ch:
        if x == 'r' or x[:-1] not in ch: continue
        p = x[:-1]; c = ch[x]; pc = ch[p]
        ip = np.array([pc['loc'][r] for r in c['rows']]); Zp = pc['z'][ip]; Zc = c['z']
        w = 1. / np.array([np.trace(pc['cov'][i]) / pc['k'] for i in ip]); w /= w.mean()
        X = np.c_[Zc, np.ones(len(Zc))]
        maps[x] = np.linalg.solve((X * w[:, None]).T @ X + 1e-9 * np.eye(X.shape[1]), (X * w[:, None]).T @ Zp)

    def compose(x, Z, Cv):
        while x != 'r':
            B = maps[x]; Z = np.c_[Z, np.ones(len(Z))] @ B; A = B[:-1]
            Cv = np.einsum('ji,njk,kl->nil', A, Cv, A); x = x[:-1]
        return Z, Cv
    def charted_ancestor(x):
        """the deepest node at or above x whose whole chain to the root has charts (deep charts may be pruned)."""
        while x != 'r':
            y = x; good = True
            while y != 'r':
                if y not in ch: good = False; break
                y = y[:-1]
            if good: return x
            x = x[:-1]
        return 'r'
    anc = {x: charted_ancestor(x) for x in set(leaf_node.values())}
    deepest = np.array([anc[leaf_node[l]] if l > 0 and l in leaf_node else 'r' for l in lab], dtype=object)
    rootloc = ch['r']['loc']; have = np.array([r in rootloc for r in range(1, N + 1)])
    kr = ch['r']['k']
    G = np.full((N, kr), np.nan); Cg = np.full((N, kr, kr), np.nan)
    for x in ch:
        sel = np.where((deepest == x) & have)[0]
        if len(sel) == 0: continue
        loc = np.array([ch[x]['loc'][r + 1] for r in sel]); Z, Cv = compose(x, ch[x]['z'][loc], ch[x]['cov'][loc]); G[sel] = Z; Cg[sel] = Cv
    Zr = np.full((N, kr), np.nan); ir = np.where(have)[0]; Zr[ir] = ch['r']['z'][[rootloc[r + 1] for r in ir]]
    eo = np.full(N, -1); eo[ir] = ch['r']['eo'][[rootloc[r + 1] for r in ir]]
    rootsplit = np.zeros(N, int)
    if ch['r']['acc'] is not None and len(ch['r']['acc']) == len(ch['r']['rows']):
        rootsplit[ch['r']['rows'] - 1] = ch['r']['acc']
    elif ch['r']['acc'] is not None and len(ch['r']['acc']) == N:
        rootsplit = ch['r']['acc'].copy()
    return dict(N=N, lab=lab, leaf_names=leaf_names, G=G, Cg=Cg, root=Zr, eo=eo, deepest=deepest, rootsplit=rootsplit, k=kr, nodes=list(ch))


def bhattacharyya(A, B):
    """overlap coefficient exp(-D_B) of Gaussians fitted to two clouds."""
    m = A.mean(0) - B.mean(0); Sa, Sb = np.cov(A.T), np.cov(B.T); S = 0.5 * (Sa + Sb) + 1e-9 * np.eye(len(m))
    db = 0.125 * m @ np.linalg.solve(S, m) + 0.5 * (np.linalg.slogdet(S)[1] - 0.5 * (np.linalg.slogdet(Sa + 1e-9 * np.eye(len(m)))[1] + np.linalg.slogdet(Sb + 1e-9 * np.eye(len(m)))[1]))
    return float(np.exp(-db))


def leaf_pairs(G, eo, lab, names, P, tag, equiv=None, allowed=None):
    """cross-half d' / AUROC and the Bhattacharyya overlap of every leaf pair (or the allowed set) in the coordinates G;
    equiv[(i, j)] = truth equivalence."""
    L = len(names); out = {}
    ok = np.isfinite(G).all(1)
    for i in range(1, L + 1):
        for j in range(i + 1, L + 1):
            if allowed is not None and (i, j) not in allowed: continue
            sel = ok & ((lab == i) | (lab == j))
            if (lab == i).sum() < 50 or (lab == j).sum() < 50: continue
            au, dp = MC.lda_crosshalf(G[sel], lab[sel] == i, eo[sel]); bc = bhattacharyya(G[ok & (lab == i)], G[ok & (lab == j)])
            out[(i, j)] = dict(dprime=dp, auroc=au, bc=bc, equiv=None if equiv is None else equiv.get((i, j)))
            P('   %-8s %-8s %s: d\' %.2f AUROC %.3f overlap %.3f%s' % (names[i - 1], names[j - 1], tag, dp, au, bc,
              '' if equiv is None else ('  [truth: %s]' % ('equivalent' if equiv.get((i, j)) else 'distinct'))))
    return out


def valleys(t, s2, eo, kmax, P, label):
    """the deconvolved 1-D density of the projections t (variances s2) fitted on each particle half: K*, the valleys (position in
    sd units of t, depth D = 1 - valley / lower neighbouring peak), the deepest valley's position per half."""
    res = {}
    sd = t.std()
    for hname, hs in (('even', eo == 0), ('odd', eo == 1)):
        K, gains, ses = MC.select_k(t[hs], s2[hs], eo[hs] * 0 + (np.arange(hs.sum()) % 2), kmax)
        f = MC.xd1(t[hs], s2[hs], K)
        lo, hi = np.quantile(t, [0.005, 0.995]); grid = np.linspace(lo, hi, 600)
        fd = MC.mix_density(f['pi'], f['xi'][:, 0], np.full(len(f['pi']), float(f['Om'][0, 0])), grid)
        lm = np.where((fd[1:-1] > fd[:-2]) & (fd[1:-1] >= fd[2:]))[0] + 1
        vs = []
        for a, b in zip(lm[:-1], lm[1:]):
            jv = a + np.argmin(fd[a:b + 1]); D = 1. - fd[jv] / max(min(fd[a], fd[b]), 1e-300)
            vs.append((float(D), float((grid[jv] - t.mean()) / sd)))
        vs.sort(reverse=True)
        res[hname] = dict(K=int(K), nmodes=int(len(lm)), valleys=vs[:3])
        P('   %s %s half: K* %d, %d modes, valleys (depth, position in sd): %s' % (label, hname, K, len(lm), ' '.join('%.2f@%+.2f' % v for v in vs[:3]) or 'none'))
    return res


def compose_cmd(a):
    T = a.tree.rstrip('/'); os.makedirs(a.out, exist_ok=True); t0 = time.time()
    rep = open(os.path.join(a.out, 'report.txt'), 'w')
    def P(s): rep.write(s + '\n'); rep.flush(); log(s)
    tr = load_tree(T); N, G, Zr, eo, lab, names = tr['N'], tr['G'], tr['root'], tr['eo'], tr['lab'], tr['leaf_names']
    kind, Tr = MC.load_truth(a.truth, N)
    P('tree %s: %d particles, %d leaves %s, charts %s, k %d' % (T, N, len(names), ' '.join(names), ' '.join(tr['nodes']), tr['k']))
    ok = np.isfinite(G).all(1) & (eo >= 0)
    summary = dict(tree=T, kind=kind, N=int(N), leaves=names)
    truth_cols = {}
    # the pairs to report: all, or (large trees) the sibling pairs, the pairs a merge file unites and the named pair
    L = len(names); allowed = None; merged_pairs = set(); named = None
    if a.pairs == 'siblings' or a.merged or a.pair:
        allowed = {(i, j) for i in range(1, L + 1) for j in range(i + 1, L + 1) if names[i - 1].split('_')[0] == names[j - 1].split('_')[0]}
        if a.pairs == 'all':
            allowed = {(i, j) for i in range(1, L + 1) for j in range(i + 1, L + 1)}
    if a.merged:
        mu = MC.read_ints(a.merged); assert len(mu) == N
        unit = {}
        for i in range(1, L + 1):
            m = (lab == i) & (mu > 0)
            if m.sum(): u, c = np.unique(mu[m], return_counts=True); unit[i] = int(u[np.argmax(c)])
        merged_pairs = {(i, j) for i in unit for j in unit if i < j and unit[i] == unit[j]}
        allowed |= merged_pairs
        P('merge file %s: %d units; %d leaf pairs united by it' % (a.merged, len(set(unit.values())), len(merged_pairs)))
    if a.pair:
        n1, n2 = a.pair.split(','); i, j = sorted((names.index(n1) + 1, names.index(n2) + 1)); named = (i, j); allowed.add(named)
    if kind == 'closed':
        closed = Tr['closedref'] > 0.5; truth_cols['closed'] = closed
    elif kind == 'qtoy':
        truth_cols = dict(conformer=Tr['conformer'].astype(int), combo=Tr['combo'].astype(int), arm=Tr['arm'], s1=Tr['s1'])
        closed = (Tr['conformer'] == 18) if 'conformer' in Tr else None
        holo = Tr['holo'] > 0.5
    # (a) minority
    if kind == 'closed' or (kind == 'qtoy' and closed is not None and closed.sum() > 0):
        P('\n== (a) minority: closed (%d of %d) vs the rest, cross-half LDA AUROC (d\')' % (closed.sum(), ok.sum()))
        for nm, Z in (('global', G), ('root', Zr)):
            au, dp = MC.lda_crosshalf(Z[ok], closed[ok], eo[ok])
            summary['closed_' + nm] = dict(auroc=au, dprime=dp)
            P('   %-7s AUROC %.3f d\' %.2f' % (nm, au, dp))
        P('   per leaf: n, closed count, enrichment; min cross-half d\' to another leaf (global / root)%s' % ('' if allowed is None else ' over the reported pairs'))
        pg = leaf_pairs(G, eo, lab, names, lambda s: None, 'global', allowed=allowed); pr = leaf_pairs(Zr, eo, lab, names, lambda s: None, 'root', allowed=allowed)
        share = closed[ok].mean()
        rows = []
        order = sorted(range(1, L + 1), key=lambda i: -closed[ok & (lab == i)].mean() if (ok & (lab == i)).sum() else 0.)
        for i in order:
            nm = names[i - 1]; m = ok & (lab == i)
            if m.sum() == 0: continue
            enr = closed[m].mean() / max(share, 1e-12)
            dg = [pg[k_]['dprime'] for k_ in pg if i in k_]; dr = [pr[k_]['dprime'] for k_ in pr if i in k_]
            rows.append(dict(leaf=nm, n=int(m.sum()), nclosed=int(closed[m].sum()), enrichment=float(enr), min_dprime_global=float(min(dg)) if dg else None, min_dprime_root=float(min(dr)) if dr else None))
            if len(rows) <= a.maxleaves:
                P('   %-12s n %7d closed %5d enrichment %5.2f | min d\' %s / %s' % (nm, m.sum(), closed[m].sum(), enr,
                  '%.2f' % min(dg) if dg else '-', '%.2f' % min(dr) if dr else '-'))
        summary['leaves'] = rows
        P('   leaf pairs in the global frame%s:' % ('' if allowed is None else ' (sibling, merged and named pairs)'))
        pg = leaf_pairs(G, eo, lab, names, P if L <= a.maxleaves else (lambda s: None), 'global', allowed=allowed)
        if L > a.maxleaves:
            ds = sorted(v['dprime'] for v in pg.values())
            P('   %d pairs: d\' quartiles %.2f / %.2f / %.2f, min %.2f max %.2f' % (len(ds), np.percentile(ds, 25), np.median(ds), np.percentile(ds, 75), ds[0], ds[-1]))
        if merged_pairs:
            dm = sorted(pg[k_]['dprime'] for k_ in pg if k_ in merged_pairs); dn = sorted(pg[k_]['dprime'] for k_ in pg if k_ not in merged_pairs)
            P('   pairs united by the merge file: d\' %s' % ' '.join('%.2f' % v for v in dm))
            P('   other reported pairs: d\' quartiles %.2f / %.2f / %.2f (%d)' % (np.percentile(dn, 25), np.median(dn), np.percentile(dn, 75), len(dn)) if dn else '   no other pairs')
            summary['merged_pairs_dprime'] = dm; summary['other_pairs_dprime'] = dn
        if named is not None and named in pg:
            v = pg[named]; rank = sum(1 for x in pg.values() if x['dprime'] > v['dprime'])
            P('   named pair %s: d\' %.2f AUROC %.3f overlap %.3f; %d of %d reported pairs are further apart; united by the merge file: %s'
              % (a.pair, v['dprime'], v['auroc'], v['bc'], rank, len(pg), named in merged_pairs))
            summary['named_pair'] = dict(pair=a.pair, **v, rank=rank, npairs=len(pg), merged=named in merged_pairs)
    # (c) leaf pairs vs truth equivalence
    if kind == 'qtoy':
        key = 'combo' if len(np.unique(truth_cols['combo'])) <= 6 else 'conformer'
        tv = truth_cols[key]
        dom = {}
        for i in range(1, len(names) + 1):
            m = ok & (lab == i)
            if m.sum() == 0: continue
            u, c = np.unique(tv[m], return_counts=True); dom[i] = (int(u[np.argmax(c)]), float(c.max() / m.sum()))
        equiv = {(i, j): dom[i][0] == dom[j][0] for i in dom for j in dom if i < j}
        P('\n== (c) leaf pairs vs truth equivalence (dominant %s per leaf: %s)' % (key, ' '.join('%s:%d(%.2f)' % (names[i - 1], d[0], d[1]) for i, d in dom.items())))
        pg = leaf_pairs(G, eo, lab, names, P, 'global', equiv)
        pr = leaf_pairs(Zr, eo, lab, names, P, 'root', equiv)
        for nm, pp_ in (('global', pg), ('root', pr)):
            de = [v['dprime'] for v in pp_.values() if v['equiv']]; dn = [v['dprime'] for v in pp_.values() if v['equiv'] is False]
            if de and dn:
                sc = np.array(de + dn); pos = np.r_[np.zeros(len(de), bool), np.ones(len(dn), bool)]
                P('   %s: d\' of equivalent pairs %s | distinct pairs %s | AUROC(d\' separates distinct from equivalent) %.3f; a merge at d\' < 1 merges %d of %d equivalent and %d of %d distinct pairs'
                  % (nm, ' '.join('%.2f' % v for v in sorted(de)), ' '.join('%.2f' % v for v in sorted(dn)), MC.auroc(sc, pos),
                     sum(v < 1. for v in de), len(de), sum(v < 1. for v in dn), len(dn)))
                summary['equiv_' + nm] = dict(equivalent=de, distinct=dn, auroc=MC.auroc(sc, pos))
    # (d) the root split's cut along the global axis
    rs = tr['rootsplit']; groups = [g for g in np.unique(rs) if g > 0]
    if len(groups) < 2:
        P('\n== (d) no root split labels (split_children.txt missing or of another length)')
    if len(groups) >= 2 and ok.sum() > 1000:
        sizes = {g: int(((rs == g) & ok).sum()) for g in groups}; g1, g2 = sorted(groups, key=lambda g: -sizes[g])[:2]
        subset = ok & (truth_cols['arm'] if 'arm' in truth_cols else True)
        if a.dsub and subset.sum() > a.dsub:
            keep = np.zeros(N, bool); keep[np.random.default_rng(0).choice(np.where(subset)[0], a.dsub, replace=False)] = True; subset &= keep
        ma, mb = subset & (rs == g1), subset & (rs == g2)
        d = G[ma].mean(0) - G[mb].mean(0); d /= np.linalg.norm(d)
        t = G[subset] @ d; s2 = np.einsum('ni,ij,nj->n', np.broadcast_to(d, (subset.sum(), len(d))), np.eye(len(d)), np.einsum('nij,j->ni', tr['Cg'][subset], d))
        s2 = np.maximum(s2, 1e-6 * t.var())
        cut = 0.5 * (np.median(G[ma] @ d) + np.median(G[mb] @ d)); sd = t.std()
        P('\n== (d) root split %d|%d (%d / %d of the %s particles) along its mean-difference direction in the global frame: the tree\'s cut at %+.2f sd; deconvolved density per half:'
          % (g1, g2, ma.sum(), mb.sum(), 'arm' if 'arm' in truth_cols else 'assigned', (cut - t.mean()) / sd))
        summary['cut'] = dict(tree_cut_sd=float((cut - t.mean()) / sd), valleys=valleys(t, s2, eo[subset], a.kmax, P, 'split axis'))
        pa = MC.pc1(G[subset]); s2p = np.full(len(pa), float(np.median(s2)))
        P('   and along PC1 of the same particles:')
        summary['cut']['pc1'] = valleys(pa, s2p, eo[subset], a.kmax, P, 'PC1')
        # the control without the composition's seams: the root chart alone (one affine frame), same particles and direction rule
        dr = Zr[ma].mean(0) - Zr[mb].mean(0); dr /= np.linalg.norm(dr)
        tr_ = Zr[subset] @ dr; cutr = 0.5 * (np.median(Zr[ma] @ dr) + np.median(Zr[mb] @ dr))
        s2r = np.full(len(tr_), float(np.median(s2)) * tr_.var() / max(t.var(), 1e-30))
        P('   root chart alone (no composition), split direction, the tree\'s cut at %+.2f sd:' % ((cutr - tr_.mean()) / tr_.std()))
        summary['cut']['root'] = dict(tree_cut_sd=float((cutr - tr_.mean()) / tr_.std()), valleys=valleys(tr_, s2r, eo[subset], a.kmax, P, 'root split axis'))
    np.savez(os.path.join(a.out, 'global.npz'), G=G, root=Zr, eo=eo, lab=lab, names=np.array(names), rootsplit=rs, **{k_: v for k_, v in truth_cols.items()})
    summary['wall'] = time.time() - t0
    json.dump(summary, open(os.path.join(a.out, 'summary.json'), 'w'), indent=1, default=float)
    P('done in %.0f s' % summary['wall']); rep.close()


def compare_cmd(a):
    A = np.load(os.path.join(a.out_a, 'global.npz'), allow_pickle=True); B = np.load(os.path.join(a.out_b, 'global.npz'), allow_pickle=True)
    assert len(A['lab']) == len(B['lab']), 'different particle counts'
    print('== (b) %s vs %s: leaf-partition NMI %.3f' % (a.name_a, a.name_b, LS.nmi(A['lab'], B['lab'])))
    def whiten(X):
        w, V = np.linalg.eigh(np.cov(X.T)); keep = w > 1e-8 * w.max()
        return X @ V[:, keep] / np.sqrt(w[keep])
    for nm in ('G', 'root'):
        Za, Zb = A[nm], B[nm]; ok = np.isfinite(Za).all(1) & np.isfinite(Zb).all(1)
        Xa = Za[ok] - Za[ok].mean(0); Xb = Zb[ok] - Zb[ok].mean(0)
        Wa, Wb = whiten(Xa), whiten(Xb)
        U, S, Vt = np.linalg.svd(Wa.T @ Wb / (ok.sum() - 1), full_matrices=False)
        cc = S
        # the first canonical pairs: per-particle agreement and extrema overlap
        rhos = []
        for c in range(min(3, len(S))):
            ua, ub = Wa @ U[:, c], Wb @ Vt[c]
            if np.corrcoef(ua, ub)[0, 1] < 0: ub = -ub
            rhos.append(MC.spearman(ua, ub))
            if c == 0:
                n5 = int(0.05 * ok.sum())
                top = len(set(np.argsort(ua)[-n5:]) & set(np.argsort(ub)[-n5:])) / n5; bot = len(set(np.argsort(ua)[:n5]) & set(np.argsort(ub)[:n5])) / n5
        print('   %-6s canonical correlations %s | canonical variates 1-3 Spearman %s; variate 1 top-5%% overlap %.2f, bottom-5%% overlap %.2f (chance 0.05)'
              % ('global' if nm == 'G' else nm, ' '.join('%.3f' % v for v in cc[:5]), ' '.join('%.3f' % r for r in rhos), top, bot))
    if 'closed' in A.files:
        cl = A['closed'] > 0.5
        for nm, X in (('A', A), ('B', B)):
            lab = X['lab']; share = cl.mean()
            enr = {i: cl[lab == i].mean() / share for i in np.unique(lab) if i > 0}
            best = max(enr, key=enr.get)
            print('   %s: most closed-enriched leaf %s (x%.2f, n %d)' % (nm, X['names'][best - 1], enr[best], (lab == best).sum()))
        la, lb = A['lab'], B['lab']
        ea = {i: cl[la == i].mean() / cl.mean() for i in np.unique(la) if i > 0}; eb = {i: cl[lb == i].mean() / cl.mean() for i in np.unique(lb) if i > 0}
        ia, ib = max(ea, key=ea.get), max(eb, key=eb.get)
        inter = ((la == ia) & (lb == ib)).sum(); print('   closed-enriched leaves overlap: %d particles of %d / %d (Jaccard %.2f)' % (inter, (la == ia).sum(), (lb == ib).sum(), inter / ((la == ia) | (lb == ib)).sum()))


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    c = sub.add_parser('compose'); c.add_argument('tree'); c.add_argument('out'); c.add_argument('--truth', default='none'); c.add_argument('--kmax', type=int, default=6)
    c.add_argument('--pairs', default='all', help='all | siblings: the leaf pairs reported (plus merged and named pairs)')
    c.add_argument('--merged', default='', help='per-particle unit labels of a merge-up (hetero_tree_merged.txt)')
    c.add_argument('--pair', default='', help='LEAF1,LEAF2: a named pair to report (the closed pair)')
    c.add_argument('--dsub', type=int, default=0, help='subsample for the density readout (d)'); c.add_argument('--maxleaves', type=int, default=12)
    b = sub.add_parser('compare'); b.add_argument('out_a'); b.add_argument('out_b'); b.add_argument('name_a'); b.add_argument('name_b')
    a = ap.parse_args()
    compose_cmd(a) if a.cmd == 'compose' else compare_cmd(a)


if __name__ == '__main__':
    main()

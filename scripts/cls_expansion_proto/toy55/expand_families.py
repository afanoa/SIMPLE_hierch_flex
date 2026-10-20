#!/usr/bin/env python3
"""Genetic expansion, step 1: split a parent klem run's states into families by spectral bisection of the state graph
(correlation of the state deviations at 60-20 A), assign every on-grid particle to the family of its argmax state and
write one subset file per family (1-based particle indices) for the child runs. Step 3 (--merge): combine the child
runs' latents into one labelling (family x surviving child state) and score coverage against the truth.
usage: expand_families.py split PARENT_DIR OUT_DIR [--nfam 2]
       expand_families.py softsplit PARENT_DIR OUT_DIR [--nfam 2 --thr 0.05 --parent_weights F]   (soft families: subsets + weights)
       expand_families.py merge OUT_DIR CHILD_DIR1 CHILD_DIR2 ... --gt GT [--minrep 0.5] [--weights w1,w2,..]"""
import sys, os, glob, re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('mode'); ap.add_argument('args', nargs='+'); ap.add_argument('--nfam', type=int, default=2); ap.add_argument('--minrep', type=float, default=0.5); ap.add_argument('--gt', default=''); ap.add_argument('--thr', type=float, default=0.05, help='softsplit: particles with a family weight below this are left out of that family'); ap.add_argument('--parent_weights', default='', help='softsplit: the parent board\'s own weights file (chained soft tree): child weights = parent weight x family responsibility'); ap.add_argument('--minsize', type=int, default=0, help='merge: states with fewer labelled particles than this are dropped before the purity readout (a 50-particle state can be pure by accident)'); ap.add_argument('--bymaxr', action='store_true', help='merge: children whose latents cover the same particles (final_all): a particle goes to the child where its max responsibility is largest (global re-assignment at the leaves)'); ap.add_argument('--weights', default='', help='merge: comma list of per-child weights files (soft tree): a particle goes to the child where weight x max responsibility is largest'); a = ap.parse_args()
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
def families(parent, nfam):
    """cell families by spectral bisection of the state graph (correlation of the cell deviations); returns fam (per cell), the similarity S and K"""
    fs = sorted(glob.glob(parent + '/cells/cell_*.mrc')) + sorted(glob.glob(parent + '/cells/state_*.mrc')); V = np.array([read_mrc(f) for f in fs]); K = len(V); n = V.shape[1]; Dv = V - V.mean(0)
    zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); band = r >= 1
    X = np.array([np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(d)))[band] for d in Dv]); X /= np.maximum(np.sqrt(np.sum(np.abs(X) ** 2, 1))[:, None], 1e-30); S = np.real(X @ np.conj(X).T)
    W = np.maximum(S, 0); np.fill_diagonal(W, 0); dg = W.sum(1) + 1e-9; Lsym = np.eye(K) - W / np.sqrt(np.outer(dg, dg)); w_, v_ = np.linalg.eigh(Lsym)
    if nfam == 2: fam = (v_[:, 1] > np.median(v_[:, 1])).astype(int)
    else:
        emb = v_[:, 1:nfam]; rng = np.random.default_rng(0); C = emb[rng.choice(K, nfam, replace=False)]
        for it in range(100): lab = ((emb[:, None, :] - C[None]) ** 2).sum(2).argmin(1); C = np.array([emb[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(nfam)])
        fam = lab
    return fam, S, K
if a.mode == 'softsplit':
    # soft families: every on-grid particle gets a weight per family = its responsibility mass on that family's cells (x the parent's own weight in a chained tree);
    # family subset = particles above --thr, weights file beside it; the children train on the union with these weights, the merge picks the best child per particle
    parent, out = a.args[:2]; os.makedirs(out, exist_ok=True); fam, S, K = families(parent, a.nfam)
    R = np.load(parent + '/klem_resp.npz'); r = R['r'].astype(float); pind = R['pinds'].astype(int); ngrid = len(R['i1']); rg = r[:, :ngrid]; gm = rg.sum(1); on = gm > 0.5; base = np.ones(len(pind))
    if a.parent_weights:
        wm = {}
        for l in open(a.parent_weights):
            if l.strip() and not l.startswith('#'): t = l.split(); wm[int(float(t[0]))] = float(t[1])
        base = np.array([wm.get(int(p_), 0.) for p_ in pind])
    for j in range(a.nfam):
        w = base * rg[:, fam[:ngrid] == j].sum(1) / np.maximum(gm, 1e-12); keep = on & (w > a.thr)
        np.savetxt('%s/family%d_subset.txt' % (out, j + 1), pind[keep], fmt='%d'); np.savetxt('%s/family%d_weights.txt' % (out, j + 1), np.c_[pind[keep], w[keep]], fmt='%d %.4f')
        print('family %d: %d cells; soft membership: %d particles above %.2f (%d above 0.5, %d above 0.9), weight sum %.0f; within-family mean state similarity %.2f' % (j + 1, (fam == j).sum(), keep.sum(), a.thr, int((w > 0.5).sum()), int((w > 0.9).sum()), w[keep].sum(), S[np.ix_(fam == j, fam == j)].mean()))
    np.savetxt(out + '/families.txt', np.c_[np.arange(1, K + 1), fam + 1], fmt='%d', header='state family')
elif a.mode == 'split':
    parent, out = a.args[:2]; os.makedirs(out, exist_ok=True)
    fs = sorted(glob.glob(parent + '/cells/cell_*.mrc')) + sorted(glob.glob(parent + '/cells/state_*.mrc')); V = np.array([read_mrc(f) for f in fs]); K = len(V); n = V.shape[1]; Dv = V - V.mean(0)
    zz, yy, xx = np.indices((n, n, n)); r = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); band = r >= 1
    X = np.array([np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(d)))[band] for d in Dv]); X /= np.maximum(np.sqrt(np.sum(np.abs(X) ** 2, 1))[:, None], 1e-30); S = np.real(X @ np.conj(X).T)
    W = np.maximum(S, 0); np.fill_diagonal(W, 0); dg = W.sum(1) + 1e-9; Lsym = np.eye(K) - W / np.sqrt(np.outer(dg, dg)); w_, v_ = np.linalg.eigh(Lsym)
    if a.nfam == 2: fam = (v_[:, 1] > np.median(v_[:, 1])).astype(int)
    else:
        emb = v_[:, 1:a.nfam]; rng = np.random.default_rng(0); C = emb[rng.choice(K, a.nfam, replace=False)]
        for it in range(100): lab = ((emb[:, None, :] - C[None]) ** 2).sum(2).argmin(1); C = np.array([emb[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(a.nfam)])
        fam = lab
    L = np.loadtxt(parent + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; on = (gm > 0.5) & (am < K); famp = np.where(am < K, fam[np.minimum(am, K - 1)], -1)
    for j in range(a.nfam):
        members = pind[on & (famp == j)]; np.savetxt('%s/family%d_subset.txt' % (out, j + 1), members, fmt='%d')
        print('family %d: %d states, %d particles; within-family mean state similarity %.2f' % (j + 1, (fam == j).sum(), len(members), S[np.ix_(fam == j, fam == j)].mean()))
    print('between-family mean state similarity %.2f' % S[np.ix_(fam == 0, fam == 1)].mean() if a.nfam == 2 else '')
    np.savetxt(out + '/families.txt', np.c_[np.arange(1, K + 1), fam + 1], fmt='%d', header='state family')
elif a.mode == 'merge':
    out = a.args[0]; children = a.args[1:]; gt = {}
    for l in open(a.gt):
        if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
    lab = {}; score = {}; nsurv = 0; ntot = 0
    for j, d in enumerate(children):
        L = np.loadtxt(d + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]
        st = np.loadtxt(d + '/klem_states.txt', comments='#'); rep = dict(zip(st[:, 0].astype(int), st[:, 2])); grid = dict(zip(st[:, 0].astype(int), st[:, 3].astype(int)))
        surv = {k for k, v in rep.items() if v >= a.minrep and grid.get(k, 0) == 1}; nsurv += len(surv); ntot += int(sum(grid.values()))
        wmap = {}
        if a.weights:
            for l in open(a.weights.split(',')[j]):
                if l.strip() and not l.startswith('#'): t = l.split(); wmap[int(float(t[0]))] = float(t[1])
        mr = L[:, 4]
        for p_, k_, g_, m_ in zip(pind, am, gm, mr):
            if g_ > 0.5 and k_ in surv:
                sc = (wmap.get(int(p_), 1.) if a.weights else 1.) * (m_ if (a.weights or a.bymaxr) else 1.)
                if int(p_) not in lab or sc > score.get(int(p_), -1.): lab[int(p_)] = (j, int(k_)); score[int(p_)] = sc
        print('child run %d (%s): %d of %d grid states survive (cross-half >= %.2f), %d particles labelled' % (j + 1, os.path.basename(d), len(surv), int(sum(grid.values())), a.minrep, sum(1 for p_ in pind if int(p_) in lab and lab[int(p_)][0] == j)))
    if not lab: print('MERGED children: no surviving states'); raise SystemExit
    if a.minsize > 0:
        cnt_ = {}
        for v in lab.values(): cnt_[v] = cnt_.get(v, 0) + 1
        small = {k for k, c in cnt_.items() if c < a.minsize}; lab = {p_: v for p_, v in lab.items() if v not in small}; print('   states below %d particles dropped: %d (%d particles left)' % (a.minsize, len(small), len(lab)))
    conf = np.array([gt.get(p_, 0) for p_ in lab]); L_ = np.array([hash(v) for v in lab.values()]); m = conf > 0; cs = sorted(set(conf[m])); us, ui = np.unique(L_[m], return_inverse=True)
    C = np.zeros((len(us), len(cs)))
    for i_, c in enumerate(cs): C[:, i_] = np.bincount(ui[conf[m] == c], minlength=len(us))
    P = C / np.maximum(C.sum(1, keepdims=True), 1); base = C.sum(0) / C.sum()
    from scipy.optimize import linear_sum_assignment
    r_, c_ = linear_sum_assignment(-P); pur = np.zeros(len(cs)); pur[c_] = P[r_, c_]; owned = len(set(P.argmax(1)))
    print('MERGED children: %d surviving states over %d particles; conformers owning a state (majority): %d/%d; one-state-per-conformer purity mean %.2f (min %.2f): %s' % (len(us), m.sum(), owned, len(cs), pur.mean(), pur.min(), ' '.join('%d:%.2f' % (c, pur[i_]) for i_, c in enumerate(cs))))
    print('   best-state purity per conformer (x base): ' + ' '.join('%d:%.2f(%.1fx)' % (c, P[:, i_].max(), P[:, i_].max() / base[i_]) for i_, c in enumerate(cs)))

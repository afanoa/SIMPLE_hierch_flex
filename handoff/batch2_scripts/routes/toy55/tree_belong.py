#!/usr/bin/env python3
"""Do the children of a data-derived tree belong to their parents? Cross-half checks on the final leaves (both halves were
split in lockstep, so a leaf name pairs the halves). (1) maps: at every ancestor split, a leaf's half-A deviation from the
parent is correlated with the half-B deviations of the alternatives there (its own branch WITHOUT the leaf, and each sibling
branch), and the reverse; the leaf belongs when its own branch wins. (2) particles: every particle is re-assigned over ALL
leaves with the OTHER half's leaf maps (no noise imprint) and with its own half's (imprinted); reported are the fractions that
change leaf / first-level branch and, with --gt, the truth purity of the tree labels against the flat re-assignment.
usage: tree_belong.py --sim A,B --poses A,B --leaves A,B --cache PREFIX [--gt A,B --truth NPZ]"""
import sys, os, argparse, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--leaves', required=True); ap.add_argument('--cache', default='')
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--cut', type=float, default=0.45); ap.add_argument('--crop', type=int, default=48)
ap.add_argument('--bandlo', type=float, default=40.); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464)
ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--skip_reassign', type=int, default=1, help='1: map test only (the particle re-assignment lives in tree_flat_estep.py)')
a = ap.parse_args()
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def setup(n):   # splitbench3.py geometry (identity conventions)
    c = a.crop; x = np.arange(n) - n // 2; rr2 = np.hypot(x[:, None], x[None, :]); R = a.mskdiam / 2 / a.smpd
    msk2 = np.clip((R + 4 - rr2) / 8., 0, 1).astype(np.float32); sl = np.r_[0:c // 2, n - c // 2:n]
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh3 = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel()
    g = np.arange(c) - c // 2; r3 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) * n / c; msk3 = np.clip((R + 4 - r3) / 8., 0, 1).astype(np.float32)
    return dict(n=n, c=c, msk2=msk2, sl=sl, fx2=fx2, fy2=fy2, kk0=kk0, sh3=sh3, msk3=msk3)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
def load_half(simdir, posef, leaff, tag):
    lab = np.loadtxt(leaff).astype(int); rows = np.where(lab > 0)[0]; names = {int(l.split()[0]): l.split()[1] for l in open(leaff + '.names') if l.strip()}
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])
    h = np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32); nx = int(h[0]); G = setup(nx); c = G['c']
    F = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r')[rows]
    pp = P[rows]; kf = np.hypot(G['fx2'], G['fy2']) / a.smpd; ang = np.arctan2(G['fy2'], G['fx2'])
    df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None, None] + (pp[:, 5] - pp[:, 6])[:, None, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None, None]))) * 1e4
    chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
    ramp = np.exp(2j * np.pi * (G['fx2'][None] * pp[:, 3, None, None] + G['fy2'][None] * pp[:, 4, None, None])).astype(np.complex64)
    Wv = (C * np.conj(ramp) * F).reshape(len(rows), -1); Wd = (C ** 2).reshape(len(rows), -1); F2 = (np.abs(F) ** 2).reshape(len(rows), -1).astype(np.float32)
    Rm = rotm(pp[:, 0], pp[:, 1], pp[:, 2]); p3 = np.einsum('qj,nij->nqi', G['kk0'], Rm); q = np.round(p3).astype(np.int32)
    ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c; idx = (q[..., 2] * c * c + q[..., 1] * c + q[..., 0]).astype(np.int64); idx[~ok] = c ** 3
    log('%s: %d particles in %d leaves' % (tag, len(rows), len(names)))
    return dict(lab=lab[rows], rows=rows, names=names, Wv=Wv, Wd=Wd, F2=F2, idx=idx, G=G)
HA = load_half(*[s.split(',')[0] for s in (a.sim, a.poses, a.leaves)], 'A'); HB = load_half(*[s.split(',')[1] for s in (a.sim, a.poses, a.leaves)], 'B')
G = HA['G']; c = G['c']; n = G['n']; sh = G['sh3']; ns = c // 2 + 1; klo = max(int(round(n * a.smpd / a.bandlo)), 3); khi = c // 2 - 1; m3 = c ** 3 + 1
nmB = {v: k for k, v in HB['names'].items()}; leaves = [HA['names'][k] for k in sorted(HA['names']) if HA['names'][k] in nmB]
log('%d leaves in both halves; bands: coarse shells 2..%d (>%.0f A), tested %d..%d (%.0f..%.1f A)' % (len(leaves), klo - 1, a.bandlo, klo, khi, n * a.smpd / klo, n * a.smpd / khi))
def acc(H, sel):
    ii = H['idx'][sel].ravel(); return np.bincount(ii, H['Wv'][sel].real.ravel(), m3)[:-1] + 1j * np.bincount(ii, H['Wv'][sel].imag.ravel(), m3)[:-1], np.bincount(ii, H['Wd'][sel].ravel(), m3)[:-1], len(sel)
for H, nm in ((HA, None), (HB, nmB)):   # per-leaf Fourier accumulators: any union of leaves (a node, a node without one leaf) is a sum
    code = {HA['names'][k]: k for k in HA['names']} if nm is None else nm; H['acc'] = {L: acc(H, np.where(H['lab'] == code[L])[0]) for L in leaves}
    H['code'] = code
nref = int(np.median([HA['acc'][L][2] for L in leaves])); den0 = 0.5 * (HA['acc'][leaves[0]][1] + HB['acc'][leaves[0]][1]) * nref / max(HA['acc'][leaves[0]][2], 1)
epsr = np.array([np.median(den0[(sh == k) & (den0 > 0)]) if ((sh == k) & (den0 > 0)).any() else 1. for k in range(ns)]) * a.eps
def vmap(H, Ls, minus=None):   # map of the union of leaves Ls (optionally without leaf 'minus'), regulariser proportional to N
    num = sum(H['acc'][L][0] for L in Ls); den = sum(H['acc'][L][1] for L in Ls); N = sum(H['acc'][L][2] for L in Ls)
    if minus: num = num - H['acc'][minus][0]; den = den - H['acc'][minus][1]; N -= H['acc'][minus][2]
    V = num / (den + epsr[sh] * max(N, 1) / nref); v = np.fft.fftshift(np.real(np.fft.ifftn(V.reshape(c, c, c)))) * G['msk3']; return np.fft.fftn(np.fft.ifftshift(v)).ravel(), N
bands = {'coarse': (sh >= 2) & (sh < klo), 'tested': (sh >= klo) & (sh <= khi)}
def corr(X, Y, b): x, y = X[bands[b]], Y[bands[b]]; return float(np.real(np.vdot(x, y)) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-30))
under = lambda nd: [L for L in leaves if L == nd or L.startswith(nd + '.')]
# ---------------- (1) map belonging at every ancestor split
log('\n(1) MAP TEST: at each split the leaf (one half) is compared with its own branch minus the leaf and with every sibling branch (other half), all minus the parent')
tally = {b: [0, 0] for b in bands}; rows_out = []
for L in leaves:
    parts = L.split('.'); line = []
    for d in range(1, len(parts)):
        par = '.'.join(parts[:d]); own = '.'.join(parts[:d + 1]); sibs = sorted(set('.'.join(x.split('.')[:d + 1]) for x in under(par)) - {own})
        if not sibs: continue
        res = {}
        for b in bands:
            sc = {}
            for X, Y in ((HA, HB), (HB, HA)):
                P_X, _ = vmap(X, under(par)); P_Y, _ = vmap(Y, under(par)); LX, _ = vmap(X, [L])
                alts = {'own': vmap(Y, under(own), minus=L if own != L else None)} if own != L else {'own(self)': vmap(Y, [L])}
                alts.update({s: vmap(Y, under(s)) for s in sibs})
                for k_, (M, N) in alts.items(): sc[k_] = sc.get(k_, 0.) + 0.5 * corr(LX - P_X, M - P_Y, b)
            win = max(sc, key=sc.get); res[b] = (win, sc)
            if own != L: tally[b][0] += win == 'own'; tally[b][1] += 1
        line.append('%s|%s: %s' % (par, 'self' if own == L else 'branch', '  '.join('%s %s' % (b, ' '.join('%s%s=%.2f' % ('*' if k_ == res[b][0] else '', k_.replace(par + '.', '.'), v) for k_, v in res[b][1].items())) for b in bands)))
    log(' %-10s n=%5d/%5d  %s' % (L, HA['acc'][L][2], HB['acc'][L][2], ' || '.join(line)))
log(' leaf closest to its own branch (minus itself) rather than a sibling, over all non-trivial ancestor splits: ' + ', '.join('%s %d/%d' % (b, v[0], v[1]) for b, v in tally.items()))
if a.skip_reassign: sys.exit(0)
# ---------------- (2) particle re-assignment over all leaves
first = lambda L: '.'.join(L.split('.')[:2])
def reassign(H, Hmap):   # E-step over all leaves with Hmap's leaf maps; Gaussian likelihood on shells 2..khi, per-shell noise from the tree labels
    Vs = [np.r_[vmap(Hmap, [L])[0], 0.].astype(np.complex64) for L in leaves]; pri = np.log(np.array([Hmap['acc'][L][2] for L in leaves], float) / sum(Hmap['acc'][L][2] for L in leaves))
    shq = np.minimum(sh[H['idx'] % (c ** 3)], ns - 1); bm = (shq >= 2) & (shq <= khi) & (H['idx'] < c ** 3); own = np.array([leaves.index(H['names'][g]) if H['names'][g] in leaves else -1 for g in H['lab']])
    R2 = np.zeros(ns); C2 = np.zeros(ns)
    for b in range(0, len(own), 4000):
        s_ = slice(b, b + 4000); P = np.stack([Vs[max(j, 0)][H['idx'][b + i]] for i, j in enumerate(own[s_])]); t = 2. * np.real(P * np.conj(H['Wv'][s_])) - H['Wd'][s_] * np.abs(P) ** 2
        u = bm[s_] & (own[s_] >= 0)[:, None]; R2 += np.bincount(shq[s_].ravel(), ((H['F2'][s_] - t) * u).ravel(), ns); C2 += np.bincount(shq[s_].ravel(), u.ravel().astype(float), ns)
    sig2 = np.maximum(R2 / np.maximum(C2, 1), 1e-12); S = np.zeros((len(own), len(leaves)), np.float32)
    for b in range(0, len(own), 4000):
        s_ = slice(b, b + 4000); ix = H['idx'][s_]; w = np.where(bm[s_], 1. / sig2[shq[s_]], 0.)
        for j in range(len(leaves)): P = Vs[j][ix]; S[s_, j] = 0.5 * ((2. * np.real(P * np.conj(H['Wv'][s_])) - H['Wd'][s_] * np.abs(P) ** 2) * w).sum(1)
    S += pri[None]; return own, S.argmax(1)
conf = {}
if a.gt:
    for H, k in ((HA, 0), (HB, 1)): conf[k] = np.loadtxt(a.gt.split(',')[k], comments='#')[:, 1].astype(int)[H['rows']]
clade = None
if a.truth:
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    z = np.load(a.truth); cl = fcluster(linkage(squareform(z['D'], checks=False), 'average'), t=a.cut, criterion='distance'); clade = dict(zip(z['ids'], cl))
def purity(tr, lab):   # share of the real particles that sit in their truth category's majority group
    m = tr > 0; t = tr[m]; g = lab[m]; tot = 0
    for v in np.unique(t): tot += np.bincount(g[t == v]).max()
    return tot / max(m.sum(), 1)
def nmi(tr, lab):
    m = tr > 0; t = np.unique(tr[m], return_inverse=True)[1]; g = np.unique(lab[m], return_inverse=True)[1]; J = np.zeros((t.max() + 1, g.max() + 1)); np.add.at(J, (t, g), 1); J /= J.sum()
    pt, pg = J.sum(1), J.sum(0); nz = J > 0; I = (J[nz] * np.log(J[nz] / np.outer(pt, pg)[nz])).sum(); Ht = -(pt[pt > 0] * np.log(pt[pt > 0])).sum(); return I / max(Ht, 1e-12)
log('\n(2) PARTICLE RE-ASSIGNMENT over all %d leaves (flat, one E-step)' % len(leaves))
fl = np.array([first(L) for L in leaves]); f_ix = {f: i for i, f in enumerate(sorted(set(fl)))}; fi = np.array([f_ix[f] for f in fl])
for H, k in ((HA, 0), (HB, 1)):
    for mode, Hm in (('other-half maps', HB if k == 0 else HA), ('own-half maps', H)):
        own, new = reassign(H, Hm); ok = own >= 0; o, nw = own[ok], new[ok]
        s = 'half %s, %s: leaf changed %.3f, first-level branch changed %.3f' % ('AB'[k], mode, np.mean(o != nw), np.mean(fi[o] != fi[nw]))
        if a.gt:
            t = conf[k][ok]; s += ' | conformer purity of first-level branches tree %.3f -> flat %.3f, of leaves %.3f -> %.3f; NMI(conformer; leaf) %.3f -> %.3f' % (purity(t, fi[o]), purity(t, fi[nw]), purity(t, o), purity(t, nw), nmi(t, o), nmi(t, nw))
            if clade is not None:
                tc = np.array([clade.get(v, 0) if v > 0 else 0 for v in t]); s += '; clade (cut %.2f) purity of first-level branches %.3f -> %.3f, NMI(clade; leaf) %.3f -> %.3f' % (a.cut, purity(tc, fi[o]), purity(tc, fi[nw]), nmi(tc, o), nmi(tc, nw))
        log(' ' + s)
        if mode == 'other-half maps':   # where each leaf's particles go
            for j, L in enumerate(leaves):
                m = o == j
                if not m.any(): continue
                dst = np.bincount(nw[m], minlength=len(leaves)) / m.sum(); top = [i for i in np.argsort(-dst)[:3] if dst[i] > 0.02]
                log('   %-10s n=%5d stays %.2f, same first-level branch %.2f; to %s' % (L, m.sum(), dst[j], dst[fi == fi[j]].sum(), ' '.join('%s:%.2f' % (leaves[i], dst[i]) for i in top)))

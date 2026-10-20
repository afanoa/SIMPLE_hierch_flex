#!/usr/bin/env python3
"""Can a flat re-assignment over ALL leaves of a data-derived tree undo early (root-level) mistakes? Every particle of one half
is scored against the OTHER half's leaf maps (no noise imprint) with the Gaussian likelihood, in two bias-free forms:
'debiased' (unmasked maps, the map-noise power sig2*den/(den+eps)^2 is subtracted from |P|^2, so a small leaf's noisier map is
not penalised) and 'equal-N' (every leaf map from the same number of particles). 'plain' (no correction) is reported to show the
size bias. Reported: fraction changing leaf / first-level branch, posterior mass outside the particle's first-level branch, and
with --gt the truth purity / NMI of the tree labels against the flat ones. usage: tree_flat_estep.py --sim A,B --poses A,B --leaves A,B --cache PREFIX [--gt A,B --truth NPZ]"""
import sys, os, argparse, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--sim', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--leaves', required=True); ap.add_argument('--cache', required=True)
ap.add_argument('--gt', default=''); ap.add_argument('--truth', default=''); ap.add_argument('--cut', type=float, default=0.45); ap.add_argument('--crop', type=int, default=48); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7)
ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--modes', default='plain,debiased,equalN'); ap.add_argument('--lim', type=float, default=0., help='E-step resolution limit in A (0 = crop Nyquist)'); a = ap.parse_args()
rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
lam = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
def load_half(simdir, posef, leaff, tag):   # splitbench3.py geometry, identity conventions
    lab = np.loadtxt(leaff).astype(int); rows = np.where(lab > 0)[0]; names = {int(l.split()[0]): l.split()[1] for l in open(leaff + '.names') if l.strip()}
    P = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(posef) if l.strip()])
    n = int(np.frombuffer(open(simdir + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
    fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    F = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r')[rows]; pp = P[rows]; kf = np.hypot(fx2, fy2) / a.smpd; ang = np.arctan2(fy2, fx2)
    df = 0.5 * ((pp[:, 5] + pp[:, 6])[:, None, None] + (pp[:, 5] - pp[:, 6])[:, None, None] * np.cos(2. * (ang[None] - np.radians(pp[:, 7])[:, None, None]))) * 1e4
    chi = np.pi * lam * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lam ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
    ramp = np.exp(2j * np.pi * (fx2[None] * pp[:, 3, None, None] + fy2[None] * pp[:, 4, None, None])).astype(np.complex64)
    Wv = (C * np.conj(ramp) * F).reshape(len(rows), -1); Wd = (C ** 2).reshape(len(rows), -1); F2 = (np.abs(F) ** 2).reshape(len(rows), -1).astype(np.float32)
    q = np.round(np.einsum('qj,nij->nqi', kk0, rotm(pp[:, 0], pp[:, 1], pp[:, 2]))).astype(np.int32); ok = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
    idx = (q[..., 2] * c * c + q[..., 1] * c + q[..., 0]).astype(np.int64); idx[~ok] = c ** 3
    log('%s: %d particles in %d leaves' % (tag, len(rows), len(names))); return dict(lab=lab[rows], rows=rows, names=names, Wv=Wv, Wd=Wd, F2=F2, idx=idx, n=n)
HA = load_half(*[s.split(',')[0] for s in (a.sim, a.poses, a.leaves)], 'A'); HB = load_half(*[s.split(',')[1] for s in (a.sim, a.poses, a.leaves)], 'B')
c = a.crop; n = HA['n']; m3 = c ** 3 + 1; ns = c // 2 + 1; khi = c // 2 - 1 if a.lim <= 0 else min(c // 2 - 1, int(round(n * a.smpd / a.lim)))
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), c // 2).ravel()
key = lambda s: [int(x) for x in s.split('.')[1:]]; nmA = {v: k for k, v in HA['names'].items()}; nmB = {v: k for k, v in HB['names'].items()}
leaves = sorted([L for L in nmA if L in nmB], key=key); first = np.array(['.'.join(L.split('.')[:2]) for L in leaves]); fix = {f: i for i, f in enumerate(sorted(set(first), key=key))}; fi = np.array([fix[f] for f in first])
for H, nm in ((HA, nmA), (HB, nmB)): H['own'] = np.array([leaves.index(H['names'][g]) if H['names'][g] in leaves else -1 for g in H['lab']])
def acc(H, sel):
    ii = H['idx'][sel].ravel(); return np.bincount(ii, H['Wv'][sel].real.ravel(), m3)[:-1] + 1j * np.bincount(ii, H['Wv'][sel].imag.ravel(), m3)[:-1], np.bincount(ii, H['Wd'][sel].ravel(), m3)[:-1]
sizes = {k: np.array([(H['own'] == j).sum() for j in range(len(leaves))]) for k, H in (('A', HA), ('B', HB))}; nref = int(np.median(np.r_[sizes['A'], sizes['B']])); nmin = int(min(sizes['A'].min(), sizes['B'].min()))
_, d0 = acc(HA, np.where(HA['own'] == 0)[0]); d0 = d0 * nref / sizes['A'][0]; epsr = np.array([np.median(d0[(sh == k) & (d0 > 0)]) if ((sh == k) & (d0 > 0)).any() else 1. for k in range(ns)]) * a.eps
log('%d leaves; equal N %d; shells 2..%d (%.1f A)' % (len(leaves), nmin, khi, n * a.smpd / khi))
def leafmaps(H, N_eq=None):   # unmasked maps (num/(den+eps), eps proportional to N) and their accumulators
    out = []
    for j in range(len(leaves)):
        sel = np.where(H['own'] == j)[0]
        if N_eq: sel = rng.choice(sel, N_eq, replace=False)
        num, den = acc(H, sel); e = epsr[sh] * len(sel) / nref; out.append((num / (den + e), den, e))
    return out
def sigma2(H, maps, own=None):   # per-pixel residual power per shell of the half's particles against their own leaf's map (from the other half)
    own = H['own'] if own is None else own; shq = np.minimum(sh[H['idx'] % (c ** 3)], ns - 1); bm = (shq >= 2) & (shq <= khi) & (H['idx'] < c ** 3) & (own >= 0)[:, None]; R2 = np.zeros(ns); C2 = np.zeros(ns)
    Vp = [np.r_[M[0], 0.] for M in maps]
    for b in range(0, len(own), 4000):
        s_ = slice(b, b + 4000); P = np.stack([Vp[max(j, 0)][H['idx'][b + i]] for i, j in enumerate(own[s_])])
        r = H['F2'][s_] - 2. * np.real(P * np.conj(H['Wv'][s_])) + H['Wd'][s_] * np.abs(P) ** 2; R2 += np.bincount(shq[s_].ravel(), (r * bm[s_]).ravel(), ns); C2 += np.bincount(shq[s_].ravel(), bm[s_].ravel().astype(float), ns)
    return np.maximum(R2 / np.maximum(C2, 1), 1e-12), shq, bm
def estep(H, maps, s2map, mode, prior, own=None):
    own = H['own'] if own is None else own; s2, shq, bm = sigma2(H, maps, own); S = np.zeros((len(own), len(maps)), np.float32)
    Vs = [np.r_[V, 0.].astype(np.complex64) for V, _, _ in maps]
    Vr = [np.r_[s2map[sh] * den / (den + e) ** 2, 0.].astype(np.float32) if mode == 'debiased' else None for V, den, e in maps]
    for b in range(0, len(own), 4000):
        s_ = slice(b, b + 4000); ix = H['idx'][s_]; w = np.where(bm[s_] | ((own[s_] < 0)[:, None] & (shq[s_] >= 2) & (shq[s_] <= khi)), 1. / s2[shq[s_]], 0.)
        for j in range(len(maps)):
            P = Vs[j][ix]; pw = np.abs(P) ** 2 - (Vr[j][ix] if Vr[j] is not None else 0.)
            S[s_, j] = 0.5 * ((2. * np.real(P * np.conj(H['Wv'][s_])) - H['Wd'][s_] * pw) * w).sum(1)
    S += np.log(prior)[None]; S -= S.max(1, keepdims=True); W = np.exp(S); W /= W.sum(1, keepdims=True); return W, s2
conf = {}; clade = None
if a.gt:
    for H, k in ((HA, 0), (HB, 1)): conf[k] = np.loadtxt(a.gt.split(',')[k], comments='#')[:, 1].astype(int)[H['rows']]
if a.truth:
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    z = np.load(a.truth); clade = dict(zip(z['ids'], fcluster(linkage(squareform(z['D'], checks=False), 'average'), t=a.cut, criterion='distance')))
def purity(tr, lab):
    m = tr > 0; t, g = tr[m], lab[m]; return sum(np.bincount(g[t == v]).max() for v in np.unique(t)) / max(m.sum(), 1)
def nmi(tr, lab):
    m = tr > 0; t = np.unique(tr[m], return_inverse=True)[1]; g = np.unique(lab[m], return_inverse=True)[1]; J = np.zeros((t.max() + 1, g.max() + 1)); np.add.at(J, (t, g), 1); J /= J.sum()
    pt, pg = J.sum(1), J.sum(0); nz = J > 0; return (J[nz] * np.log(J[nz] / np.outer(pt, pg)[nz])).sum() / max(-(pt[pt > 0] * np.log(pt[pt > 0])).sum(), 1e-12)
if a.gt and 'truth' in a.modes:
    log('\n== truth ceiling: the same E-step with maps of the TRUE groups (%s + junk) from the other half' % ('clades at cut %.2f' % a.cut if clade is not None else 'conformers'))
    tg = {}
    for H, k in ((HA, 0), (HB, 1)):
        t_ = conf[k]; tg[k] = np.array([0 if v == 0 else (clade.get(v, -1) if clade is not None else v) for v in t_])
    ids = sorted(set(tg[0]) | set(tg[1])); ids = [v for v in ids if v >= 0 and (tg[0] == v).sum() >= 50 and (tg[1] == v).sum() >= 50]
    for H, k in ((HA, 0), (HB, 1)):
        Hm, km = (HB, 1) if k == 0 else (HA, 0); mp_ = []
        for v in ids:
            sel = np.where(tg[km] == v)[0]; num, den = acc(Hm, sel); e = epsr[sh] * len(sel) / nref; mp_.append((num / (den + e), den, e))
        ownt = np.array([ids.index(v) if v in ids else -1 for v in tg[k]]); prior = np.array([(tg[km] == v).sum() for v in ids], float); prior /= prior.sum()
        W, _ = estep(H, mp_, None, 'plain', prior, ownt); asg = W.argmax(1); ok = ownt >= 0
        log(' half %s: truth-group accuracy %.3f (chance by the largest group %.3f); NMI(truth group; assigned) %.3f; per group recall %s' % ('AB'[k], np.mean(asg[ok] == ownt[ok]), np.bincount(ownt[ok]).max() / ok.sum(), nmi(ownt[ok] + 1, asg[ok]), ' '.join('%s:%.2f' % ('junk' if v == 0 else v, np.mean(asg[ownt == i] == i)) for i, v in enumerate(ids))))
s2half = {}
for H, k in ((HA, 0), (HB, 1)):   # each half's own image-noise power, for the other half's map-noise term
    s2half[k] = sigma2(H, leafmaps(HB if k == 0 else HA))[0]
for mode in [m_ for m_ in a.modes.replace(':', ',').split(',') if m_ != 'truth']:
    log('\n== %s' % mode)
    for H, k in ((HA, 0), (HB, 1)):
        Hm = HB if k == 0 else HA; maps = leafmaps(Hm, nmin if mode == 'equalN' else None); prior = sizes['AB'[1 - k]] / sizes['AB'[1 - k]].sum()
        if mode == 'filtered':   # consensus + Wiener-weighted leaf deviation: weight per shell = cross-half FSC of that leaf's deviation (reference noise suppressed where it does not reproduce)
            mine = leafmaps(H); cons = []
            for M in (maps, mine):
                num = sum(V * (d + e) for V, d, e in M); den = sum(d for _, d, _ in M); cons.append(num / (den + epsr[sh] * sum(sizes['A']) / nref))
            fmaps = []
            for j in range(len(leaves)):
                D, Do = maps[j][0] - cons[0], mine[j][0] - cons[1]
                f = np.array([np.real(np.vdot(D[sh == q], Do[sh == q])) / max(np.sqrt(np.vdot(D[sh == q], D[sh == q]).real * np.vdot(Do[sh == q], Do[sh == q]).real), 1e-30) for q in range(ns)])
                f = np.clip(np.convolve(np.pad(f, 1, mode='edge'), np.ones(3) / 3., mode='valid'), 0., 1.); fmaps.append((cons[0] + f[sh] * D, maps[j][1], maps[j][2]))
                if k == 0 and j < 30: log('   %-10s deviation FSC (shells 2,4,6,8,12,16,20): %s' % (leaves[j], ' '.join('%.2f' % f[q] for q in (2, 4, 6, 8, 12, 16, 20) if q < ns)))
            maps = fmaps
        W, _ = estep(H, maps, s2half[1 - k], 'plain' if mode == 'filtered' else mode, prior); own = H['own']; ok = own >= 0; o = own[ok]; nw = W.argmax(1)[ok]; Wf = np.stack([W[:, fi == f].sum(1) for f in range(fi.max() + 1)], 1)[ok]
        s = ' half %s: leaf changed %.3f, first-level branch changed %.3f, posterior mass outside the own first-level branch %.3f' % ('AB'[k], np.mean(o != nw), np.mean(fi[o] != fi[nw]), 1. - Wf[np.arange(len(o)), fi[o]].mean())
        if a.gt:
            t = conf[k][ok]; s += '\n   conformer purity: first-level branches tree %.3f -> flat %.3f, leaves %.3f -> %.3f; NMI(conformer; leaf) %.3f -> %.3f, NMI(conformer; first-level branch) %.3f -> %.3f' % (purity(t, fi[o]), purity(t, fi[nw]), purity(t, o), purity(t, nw), nmi(t, o), nmi(t, nw), nmi(t, fi[o]), nmi(t, fi[nw]))
            if clade is not None:
                tc = np.array([clade.get(v, 0) if v > 0 else 0 for v in t]); s += '\n   clade (cut %.2f) purity: first-level branches %.3f -> %.3f; NMI(clade; leaf) %.3f -> %.3f, NMI(clade; first-level branch) %.3f -> %.3f' % (a.cut, purity(tc, fi[o]), purity(tc, fi[nw]), nmi(tc, o), nmi(tc, nw), nmi(tc, fi[o]), nmi(tc, fi[nw]))
            jt = t == 0; Jl = [j_ for j_ in range(len(leaves)) if (o == j_).any() and np.mean(t[o == j_] == 0) > 0.5]
            if jt.any() and Jl: s += '\n   junk particles inside the junk-majority leaves %s: tree %.3f -> flat %.3f' % (','.join(leaves[j_] for j_ in Jl), np.isin(o[jt], Jl).mean(), np.isin(nw[jt], Jl).mean())
        log(s)
        for j, L in enumerate(leaves):
            m = o == j
            if not m.any(): continue
            dst = np.bincount(nw[m], minlength=len(leaves)) / m.sum(); top = [i for i in np.argsort(-dst)[:3] if dst[i] > 0.02]
            log('   %-10s n=%5d stays %.2f, same first-level branch %.2f; to %s' % (L, m.sum(), dst[j], dst[fi == fi[j]].sum(), ' '.join('%s:%.2f' % (leaves[i], dst[i]) for i in top)))

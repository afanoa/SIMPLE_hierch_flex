#!/usr/bin/env python3
"""land_real.py: TRUTH-FREE occupancy along one tree node's leading axis on REAL data (land_e1b.py's validated steps, no truth).
(1) START: leaf maps of the node (each half stack, crop 48) -> cross-half PCA of the leaf deviations -> every particle starts at its leaf's PC1
    score. (2) COORDINATE: NIT passes of the per-particle least-squares position lam along the end-to-end direction of the OTHER half
    (pixel sets A and B separately; split z calibrates the noise). (3) HELD-OUT WINDOWS: windows on lam_A, maps from pixel set B (and the swap)
    -> NNLS composition on each PHYSICAL basis (a rigid swing: --swing DIR with truth_swing.txt; deposited conformers: --confs DIR),
    mean lam of the test pixels -> monotone calibration lam(t) on the swing, window noise scale c_win; map-level occupancy on each basis.
    (4) Labels for SIMPLE reconstructions: NWREC equal-count windows of lam per master row (0 = outside the node).
usage: land_real.py --tree R --node r1|all --dump POSES --sims SA,SB --mr MRA,MRB --cache CP --swing DIR --confs DIR[,DIR] --out OUT"""
import os, sys, glob, argparse, json, numpy as np
from scipy.optimize import nnls
ap = argparse.ArgumentParser()
for k in ('tree', 'node', 'dump', 'sims', 'mr', 'cache', 'swing', 'confs', 'out'): ap.add_argument('--' + k, required=True)
ap.add_argument('--nit', type=int, default=2); ap.add_argument('--nw', type=int, default=10); ap.add_argument('--nwrec', type=int, default=12); ap.add_argument('--qend', type=float, default=15.)
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=232.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=1); ap.add_argument('--eps', type=float, default=0.05)
ap.add_argument('--nview', type=int, default=40); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--bstep', type=int, default=1)
ap.add_argument('--init', default='pca', help='pca, or an arm_readout table of the leaf maps in leaves.txt order: start at the leaf\'s closed-site occupancy')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n_ = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n_ ** 3).reshape(n_, n_, n_).astype(np.float64)
# ---------------- node rows, leaf labels, poses
L = [l.split() for l in open(a.tree + '/leaves.txt') if l.strip()]; LAB = np.loadtxt(a.tree + '/final_labels.txt', dtype=int); NM = len(LAB)
leaves = [j for j, l in enumerate(L) if a.node == 'all' or l[0].startswith(a.node + '_') or (l[0][:len(a.node)] == a.node and l[0][len(a.node):len(a.node) + 1].isdigit())]
junk = [j for j, l in enumerate(L) if len(l) > 7 and l[7] == 'junk']; leaves = [j for j in leaves if j not in junk]
inode = np.isin(LAB - 1, leaves) & (LAB > 0); rows_all = np.where(inode)[0]
log('node %s: %d leaves (%s), %d particles of %d master rows' % (a.node, len(leaves), ' '.join(L[j][0] for j in leaves), len(rows_all), NM))
need = np.zeros(NM, bool); need[rows_all] = True; P_all = np.zeros((NM, 8)); EOd = np.full(NM, -1)
with open(a.dump) as f:
    for r, l in enumerate(f):
        if r < NM and need[r]:
            d = dict(t.split('=', 1) for t in l.split() if '=' in t); P_all[r] = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]; EOd[r] = int(float(d.get('eo', -1)))
log('poses read from %s' % a.dump)
# ---------------- particle Fourier data per half stack
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sims.split(','); MRS = a.mr.split(','); n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]; psh = np.round(r2d[pcols]).astype(int)
ppar = np.mod(np.round(kk0[pcols, 0]).astype(int) + np.round(kk0[pcols, 1]).astype(int), 2) == 1   # True = pixel set B
klp = int(round(n * a.smpd / a.lp)); pin = (psh >= 2) & (psh <= klp)
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
ns_ = c // 2 + 1; khi = c // 2 - 1; KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
FB = (sh >= 2) & (sh <= min(klp, khi)); COLS = {'all': np.arange(len(pcols)), 'A': np.where(~ppar)[0], 'B': np.where(ppar)[0]}
def load_half(h):
    mr = np.loadtxt(MRS[h]).astype(int) - 1; crow = np.where(inode[mr])[0]; mrow = mr[crow]; Fc = np.load('%s_%s_c%d.npy' % (a.cache, 'AB'[h], c), mmap_mode='r'); pp = P_all[mrow]; nn = len(crow)
    Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32); pw = np.zeros((2, ns_)); pc = np.zeros((2, ns_))
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]; F2s = np.zeros(nn)
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); p = pp[s_]; F = np.asarray(Fc[crow[s_]]).reshape(len(p), -1)[:, pcols]
        df = 0.5 * ((p[:, 5] + p[:, 6])[:, None] + (p[:, 5] - p[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(p[:, 7])[:, None]))) * 1e4
        chi = np.pi * lamb * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lamb ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[pcols][None] * p[:, 3, None] + fy2.ravel()[pcols][None] * p[:, 4, None])).astype(np.complex64)
        Wv[s_] = C * np.conj(ramp) * F; Wd[s_] = C ** 2; F2 = np.abs(F) ** 2
        for q_ in (0, 1): m_ = ppar == bool(q_); pw[q_] += np.bincount(psh[m_], F2[:, m_].sum(0), ns_); pc[q_] += np.bincount(psh[m_], minlength=ns_) * len(p)
        q = np.round(np.einsum('qj,nij->nqi', kk0[pcols], rotm(p[:, 0], p[:, 1], p[:, 2]))).astype(np.int32); okq = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
        ii = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ii[~okq] = c ** 3; idx[s_] = ii
    wpx = 1. / np.maximum(pw / np.maximum(pc, 1), 1e-30); wpix = np.where(ppar, wpx[1][psh], wpx[0][psh]).astype(np.float32)
    ppn = np.zeros(nn, np.float32)   # per-particle noise level: own image power over the half's per-shell average in the lam band
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); ppn[s_] = (np.abs(np.asarray(Fc[crow[s_]]).reshape(len(crow[s_]), -1)[:, pcols][:, pin]) ** 2 * wpix[pin][None]).mean(1)
    ph, th = np.radians(pp[:, 0]), np.radians(pp[:, 1]); v = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    i = np.arange(a.nview) + 0.5; z = 1 - 2 * i / a.nview; rr = np.sqrt(np.maximum(1 - z * z, 0)); gg = np.pi * (3 - np.sqrt(5)) * i
    vb = np.argmax(v @ np.c_[rr * np.cos(gg), rr * np.sin(gg), z].T, 1); dfm = 0.5 * (pp[:, 5] + pp[:, 6]); dq = np.searchsorted(np.quantile(dfm, [1 / 3., 2 / 3.]), dfm)
    return dict(Wv=Wv, Wd=Wd, idx=idx, mrow=mrow, st=vb * 3 + dq, n=nn, wpix=wpix, ppn=ppn, leaf=LAB[mrow] - 1)
HS = [load_half(h) for h in (0, 1)]; N = HS[0]['n'] + HS[1]['n']
log('map data: halves %d / %d particles, crop %d, %d pixels per image, lam band shells 2..%d (%.1f A), ctfabs %d' % (HS[0]['n'], HS[1]['n'], c, len(pcols), klp, n * a.smpd / klp, a.ctfabs))
neq0 = 5000
def rec_raw(H, s_, cs='all'):
    cc = COLS[cs]; ii = H['idx'][np.ix_(s_, cc)].ravel(); m = c ** 3 + 1; Wv_ = H['Wv'][np.ix_(s_, cc)]
    return np.bincount(ii, Wv_.real.ravel(), m)[:-1] + 1j * np.bincount(ii, Wv_.imag.ravel(), m)[:-1], np.bincount(ii, H['Wd'][np.ix_(s_, cc)].ravel(), m)[:-1]
def shell_eps(H):
    den = rec_raw(H, np.arange(H['n']))[1] * (neq0 / H['n']); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns_)]) * a.eps
EPS = 0.5 * (shell_eps(HS[0]) + shell_eps(HS[1]))
def rec(H, s_, cs='all'):
    num, den = rec_raw(H, s_, cs); f_ = len(COLS[cs]) / len(pcols); V_ = num / (den + EPS[sh] * f_ * max(len(s_), 1) / neq0)
    v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
def dock(Mp, Ct, m_=FB):
    a_ = Mp[m_] * np.conj(Ct[m_]); k_ = KV[m_]; t = np.zeros(3)
    for _ in range(8):
        z = a_ * np.exp(-2j * np.pi * (k_ @ t)); g_ = 2 * np.pi * (k_ * z.imag[:, None]).sum(0); H_ = -(2 * np.pi) ** 2 * (k_.T * z.real) @ k_
        try: st_ = np.linalg.solve(H_, g_)
        except np.linalg.LinAlgError: break
        t = t - st_
        if np.abs(st_).max() < 1e-4: break
    if not np.all(np.isfinite(t)) or np.abs(t).max() > 3.: return Mp
    return Mp * np.exp(-2j * np.pi * (KV @ t))
CONS = [rec(H, np.arange(H['n'])) for H in HS]
# ---------------- start: cross-half PCA of the node's leaf maps
LM = [[rec(H, np.where(H['leaf'] == j)[0]) for j in leaves] for H in HS]; cnt = np.array([[(H['leaf'] == j).sum() for j in leaves] for H in HS]).sum(0)
Dv = [np.array([dock(m, CONS[h]) - CONS[h] for m in LM[h]])[:, FB] for h in (0, 1)]
Cx = np.real(Dv[0] @ np.conj(Dv[1]).T); Cx = 0.5 * (Cx + Cx.T); w_ = np.sqrt(cnt / cnt.sum()); ev, EV = np.linalg.eigh(Cx * w_[:, None] * w_[None]); pc1 = EV[:, -1] / w_
log('leaf-map PCA (cross-half): leading shares %s; PC1 leaf scores %s' % (' '.join('%.2f' % x for x in ev[::-1][:4] / max(ev[ev > 0].sum(), 1e-30)),
    ' '.join('%s:%+.2f' % (L[j][0], s) for j, s in zip(leaves, pc1))))
z = np.zeros(N); off = [0, HS[0]['n']]; sc0 = pc1
if a.init != 'pca':   # start along the physical readout: each leaf's closed-site occupancy
    occ = np.array([float(l.split()[1]) for l in open(a.init) if l.strip() and not l.startswith('#')]); sc0 = occ[np.array(leaves)]
    log('start: leaf closed-site occupancy from %s (%.2f..%.2f)' % (a.init, sc0.min(), sc0.max()))
for h, H in enumerate(HS): z[off[h]:off[h] + H['n']] = sc0[np.searchsorted(leaves, H['leaf'])] + 1e-6 * rng.normal(size=H['n'])
z = (z - z.mean()) / z.std()
def beta_of(H, M0):
    bn = np.zeros(ns_); bd = np.zeros(ns_); M0x = np.r_[M0, 0]
    for b in range(0, H['n'], 20000):
        s_ = slice(b, b + 20000); P0 = M0x[H['idx'][s_]]
        bn += np.bincount(psh, np.real(P0 * np.conj(H['Wv'][s_])).sum(0), ns_); bd += np.bincount(psh, (H['Wd'][s_] * np.abs(P0) ** 2).sum(0), ns_)
    return np.clip(bn / np.maximum(bd, 1e-30), 0., None)
BETA = [beta_of(HS[h], CONS[1 - h]) for h in (0, 1)]
def lam_of(H, h, D, M0):
    Dx = np.r_[D, 0]; M0x = np.r_[M0, 0]; b = BETA[h][psh][None]; out = np.zeros((H['n'], 4))
    for b0 in range(0, H['n'], 20000):
        s_ = slice(b0, b0 + 20000); PD = b * Dx[H['idx'][s_]]; P0 = b * M0x[H['idx'][s_]]; w = H['wpix'][None] / np.maximum(H['ppn'][s_, None], 1e-6) * pin[None]
        r = H['Wv'][s_] - H['Wd'][s_] * P0; num = w * np.real(np.conj(PD) * r); den = w * H['Wd'][s_] * np.abs(PD) ** 2
        nA, dA, nB, dB = num[:, ~ppar].sum(1), den[:, ~ppar].sum(1), num[:, ppar].sum(1), den[:, ppar].sum(1); out[s_] = np.c_[nA / np.maximum(dA, 1e-30), dA, nB / np.maximum(dB, 1e-30), dB]
    return out
def hsel(h, mask): return np.where(mask[off[h]:off[h] + HS[h]['n']])[0]
for it in range(1, a.nit + 1):
    qlo, qhi = np.percentile(z, [a.qend, 100 - a.qend]); LH = []
    for h, H in enumerate(HS): LH.append((dock(rec(H, hsel(h, z <= qlo)), CONS[h]), dock(rec(H, hsel(h, z >= qhi)), CONS[h])))
    SS = np.zeros((N, 4))
    for h, H in enumerate(HS): SS[off[h]:off[h] + H['n']] = lam_of(H, h, LH[1 - h][1] - LH[1 - h][0], CONS[1 - h])
    I = SS[:, 1] + SS[:, 3]; lam = (SS[:, 0] * SS[:, 1] + SS[:, 2] * SS[:, 3]) / np.maximum(I, 1e-30); zsp = (SS[:, 0] - SS[:, 2]) / np.sqrt(1. / np.maximum(SS[:, 1], 1e-30) + 1. / np.maximum(SS[:, 3], 1e-30))
    log('  pass %d: split z var %.3f, Fisher sd median %.3f, lam sd %.3f (end-to-end units); corr with the start %.3f' % (it, zsp.var(), np.median(1. / np.sqrt(np.maximum(I, 1e-30))), lam.std(), np.corrcoef(lam, z)[0, 1]))
    z = (lam - lam.mean()) / lam.std()
lamA, IA, lamB, IB = SS[:, 0].copy(), SS[:, 1].copy(), SS[:, 2].copy(), SS[:, 3].copy()
# ---------------- physical bases
def crop_four(v):
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; Fv = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(Fv[np.ix_(sb, sb, sb)]))) * msk3
    return np.fft.fftn(np.fft.ifftshift(vc)).ravel()
tr = [l.split() for l in open(a.swing + '/truth_swing.txt') if l.strip() and not l.startswith('#')][::a.bstep]; tips = np.array([float(l[3]) for l in tr])
SWB = np.array([crop_four(rd(glob.glob('%s/conf%s_*.mrc' % (a.swing, l[1]))[0])) for l in tr])
CF = sorted(sum([glob.glob(d + '/conf*_m160.mrc') for d in a.confs.split(',')], [])); CFB = np.array([crop_four(rd(f)) for f in CF]); CFN = [os.path.basename(f).split('_')[0] for f in CF]
def unmix(Mw, BAS):
    BM = BAS.mean(0); Mw = dock(Mw, BM); m_ = FB; NB = len(BAS); q = np.full(NB, 1. / NB)
    for _ in range(6):
        Bq = q @ BAS; aa = np.bincount(sh[m_], np.real(Mw[m_] * np.conj(Bq[m_])), ns_) / np.maximum(np.bincount(sh[m_], np.abs(Bq[m_]) ** 2, ns_), 1e-30)
        Amat = (BAS[:, m_] * aa[sh[m_]][None]).T; q, rn = nnls(np.r_[Amat.real, Amat.imag], np.r_[Mw[m_].real, Mw[m_].imag])
        if q.sum() <= 0: q = np.full(NB, 1. / NB)
    return q / q.sum(), rn / np.linalg.norm(np.r_[Mw[m_].real, Mw[m_].imag])
# ---------------- held-out windows
FOLDS = (('A', 'B'), ('B', 'A')); LS = {'A': (lamA, IA), 'B': (lamB, IB)}; out = {}
for f_, (fs, ft) in enumerate(FOLDS):
    lsel, ltest, Itest = LS[fs][0], LS[ft][0], LS[ft][1]
    qs = np.percentile(lsel, np.linspace(0, 100, a.nw + 1)); wl = np.clip(np.searchsorted(qs, lsel, side='right') - 1, 0, a.nw - 1)
    for w in range(a.nw):
        m_ = wl == w; rec_ = {}
        for h, H in enumerate(HS):
            Mw = rec(H, hsel(h, m_), ft); qs_, rs_ = unmix(Mw, SWB); qc_, rc_ = unmix(Mw, CFB); rec_[h] = (qs_, rs_, qc_, rc_)
        out[(f_, w)] = dict(n=int(m_.sum()), lam=float(ltest[m_].mean()), var=float(ltest[m_].var()), iinv=float(np.mean(1. / np.maximum(Itest[m_], 1e-30))),
                            qs=[rec_[h][0].tolist() for h in (0, 1)], rs=[rec_[h][1] for h in (0, 1)], qc=[rec_[h][2].tolist() for h in (0, 1)], rc=[rec_[h][3] for h in (0, 1)])
        o = out[(f_, w)]; qm = 0.5 * (np.array(o['qs'][0]) + np.array(o['qs'][1])); qcm = 0.5 * (np.array(o['qc'][0]) + np.array(o['qc'][1])); tc = np.argsort(-qcm)[:3]
        log('  fold %s->%s window %2d: %6d ptcls, lam(test) %+.3f | swing mean %.1f A (halves %.1f / %.1f, resid %.3f) | conformers %s (resid %.3f)' % (fs, ft, w + 1, o['n'], o['lam'],
            qm @ tips, np.array(o['qs'][0]) @ tips, np.array(o['qs'][1]) @ tips, np.mean(o['rs']), ' '.join('%s:%.2f' % (CFN[k], qcm[k]) for k in tc), np.mean(o['rc'])))
# orientation (truth-free): lam increases with the swing travel
LWall = np.array([out[k]['lam'] for k in sorted(out)]); TWall = np.array([0.5 * (np.array(out[k]['qs'][0]) + np.array(out[k]['qs'][1])) @ tips for k in sorted(out)])
sgn = np.sign(np.corrcoef(LWall, TWall)[0, 1]); log('orientation: corr(window lam, swing travel) %+.3f' % np.corrcoef(LWall, TWall)[0, 1])
if sgn < 0: lam, lamA, lamB, zsp = -lam, -lamA, -lamB, -zsp
for k in out: out[k]['lam'] *= sgn
# ---------------- labels for SIMPLE reconstructions + save
qs = np.percentile(lam, np.linspace(0, 100, a.nwrec + 1)); wr = np.clip(np.searchsorted(qs, lam, side='right') - 1, 0, a.nwrec - 1) + 1
labm = np.zeros(NM, int); mrow_all = np.concatenate([H['mrow'] for H in HS]); labm[mrow_all] = wr; np.savetxt(a.out + '/windows_rec.txt', labm, fmt='%d')
np.savez(a.out + '/real_particles.npz', mrow=mrow_all, half=np.r_[np.zeros(HS[0]['n'], int), np.ones(HS[1]['n'], int)], lam=lam, I=I, lamA=lamA, IA=IA, lamB=lamB, IB=IB, zsplit=zsp,
         st=np.concatenate([H['st'] for H in HS]), leaf=np.concatenate([H['leaf'] for H in HS]), tips=tips, conf_names=np.array(CFN))
json.dump(dict(node=a.node, n=int(N), leaves=[L[j][0] for j in leaves], pc1=pc1.tolist(), windows={'%d_%d' % k: v for k, v in out.items()}, tips=tips.tolist(), conf_names=CFN,
               c_split=float(zsp.var())), open(a.out + '/real.json', 'w'), indent=1)
log('wrote %s/real.json, real_particles.npz, windows_rec.txt (%d windows for reconstruction)' % (a.out, a.nwrec)); log('DONE %s' % a.out)

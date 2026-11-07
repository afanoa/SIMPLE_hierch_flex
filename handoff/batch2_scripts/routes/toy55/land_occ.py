#!/usr/bin/env python3
"""land_occ.py: occupancy landscape in PHYSICAL units = a linear map readout (closed-site and open-site occupancy, arm_readout.py's
definition evaluated on crop-48 maps). Held-out windows (windows on lam of pixel set S, maps from the other set, both half stacks) give each
window's mean occupancy; window means calibrate the per-particle coordinate (b2 = d mean lam_test / d mean occ), the shrinkage b1 =
d mean occ / d mean lam_S gives its reliability b1 b2 and the noise scale c_reg; the smooth cross-half deconvolution (land_deconv) gives
p(occ) with bootstrap bands. Synthetic runs (--kind e1b) are scored against occ(t_i) of the truth maps read out the same way.
land_joint.py (origin of the data code): JOINT map + particle fit of an occupancy landscape p(t) on a physical basis grid (t_k, maps V_k), the coordinate's calibration
lam(t) (monotone, knots) and its noise scale c. Per fold (windows on lam of pixel set S, maps from the other pixel set) and half stack:
  particle term: lam_S,i ~ sum_k p_k N(lam(t_k), c / I_S,i)
  map term:      window map M_w ~ a(shell) * sum_k qhat_w(k) V_k, qhat_w(k) = mean over the window's particles of their model posterior
                 P(k | lam_S,i); per-shell noise from the two half stacks' difference; a(shell) profiled out
  prior:         lam_pen ||D2 G||^2, G = -ln p
The window maps measure how much the window compositions overlap, i.e. the noise the deconvolution needs: they constrain c and lam(t),
which the particle term alone cannot separate from the width of p. Input: a land_e1b.py (synthetic, --gt) or land_real.py run directory.
usage: land_joint.py --run DIR --kind e1b|real [data args as the run] --basis LIST|--swing DIR --out OUT [--nw 12 --lam_pen 1e-3,1e-2 --nboot 10]"""
import os, sys, glob, argparse, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import land_deconv as LD
ap = argparse.ArgumentParser()
for k in ('run', 'kind', 'sim', 'mr', 'cache', 'poses', 'out'): ap.add_argument('--' + k, required=True)
ap.add_argument('--cons', required=True); ap.add_argument('--closed', required=True); ap.add_argument('--truthmaps', default=''); ap.add_argument('--grid', default='-0.4:1.4:37'); ap.add_argument('--gt', default=''); ap.add_argument('--codes', default='swing:200:230')
ap.add_argument('--nw', type=int, default=12); ap.add_argument('--nboot', type=int, default=20)
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05)
ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--tmax', type=float, default=-1., help='real: use basis steps with travel <= TMAX')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n_ = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n_ ** 3).reshape(n_, n_, n_).astype(np.float64)
# ---------------- run data: particles (master rows, half, lam per pixel set, Fisher), truth if synthetic
if a.kind == 'e1b':
    Z = np.load(a.run + '/e1b_particles.npz'); mrow = Z['rows']; half = Z['eo'].astype(int); tt = Z['t']
else:
    Z = np.load(a.run + '/real_particles.npz'); mrow = Z['mrow']; half = Z['half'].astype(int); tt = np.full(len(mrow), -1)
lamA, IA, lamB, IB = Z['lamA'], Z['IA'], Z['lamB'], Z['IB']; N = len(mrow)
# ---------------- Fourier data (crop c), halves = eo (one stack) or the two half stacks
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]; psh = np.round(r2d[pcols]).astype(int)
ppar = np.mod(np.round(kk0[pcols, 0]).astype(int) + np.round(kk0[pcols, 1]).astype(int), 2) == 1
klp = int(round(n * a.smpd / a.lp))
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
ns_ = c // 2 + 1; khi = c // 2 - 1; KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
FB = (sh >= 2) & (sh <= min(klp, khi)); COLS = {'A': np.where(~ppar)[0], 'B': np.where(ppar)[0]}
if a.kind == 'e1b':
    NMp = sum(1 for _ in open(a.poses)); need = np.zeros(NMp, bool); need[mrow] = True; P_all = np.zeros((NMp, 8))
    for r, l in enumerate(open(a.poses)):
        if need[r]: d = kv(l); P_all[r] = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]
    HR = [(np.where(half == h)[0], None) for h in (0, 1)]   # (particle indices, cache rows = master rows)
else:
    NMp = int(mrow.max()) + 1; need = np.zeros(NMp, bool); need[mrow] = True; P_all = np.zeros((NMp, 8))
    for r, l in enumerate(open(a.poses)):
        if r < NMp and need[r]: d = kv(l); P_all[r] = [float(d.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')]
    HR = []
    for h in (0, 1):
        mr = np.loadtxt(a.mr.split(',')[h]).astype(int) - 1; pos = {m: i for i, m in enumerate(mr)}; pi_ = np.where(half == h)[0]; HR.append((pi_, np.array([pos[m] for m in mrow[pi_]])))
def load_half(h):
    pi_, crow = HR[h]; mr_ = mrow[pi_]; crow = mr_ if crow is None else crow; Fc = np.load('%s_%s_c%d.npy' % (a.cache, 'A' if a.kind == 'e1b' else 'AB'[h], c), mmap_mode='r')
    pp = P_all[mr_]; nn = len(pi_); Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32)
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]; o = np.argsort(crow); inv = np.empty_like(o); inv[o] = np.arange(nn)
    for b in range(0, nn, 20000):
        s_ = slice(b, b + 20000); p = pp[s_]; F = np.asarray(Fc[crow[s_]]).reshape(len(p), -1)[:, pcols]
        df = 0.5 * ((p[:, 5] + p[:, 6])[:, None] + (p[:, 5] - p[:, 6])[:, None] * np.cos(2. * (ang[None] - np.radians(p[:, 7])[:, None]))) * 1e4
        chi = np.pi * lamb * df * kf[None] ** 2 - 0.5 * np.pi * a.cs * 1e7 * lamb ** 3 * kf[None] ** 4; C = (np.sqrt(1. - a.fraca ** 2) * np.sin(chi) + a.fraca * np.cos(chi)).astype(np.float32)
        if a.ctfabs: C = np.abs(C)
        ramp = np.exp(2j * np.pi * (fx2.ravel()[pcols][None] * p[:, 3, None] + fy2.ravel()[pcols][None] * p[:, 4, None])).astype(np.complex64)
        Wv[s_] = C * np.conj(ramp) * F; Wd[s_] = C ** 2
        q = np.round(np.einsum('qj,nij->nqi', kk0[pcols], rotm(p[:, 0], p[:, 1], p[:, 2]))).astype(np.int32); okq = (np.abs(q) <= c // 2 - 1).all(2); q = q % c
        ii = q[..., 2] * c * c + q[..., 1] * c + q[..., 0]; ii[~okq] = c ** 3; idx[s_] = ii
    return dict(Wv=Wv, Wd=Wd, idx=idx, pi=pi_, n=nn)
HS = [load_half(h) for h in (0, 1)]; log('%s: %d particles (halves %d / %d)' % (a.run, N, HS[0]['n'], HS[1]['n']))
neq0 = 5000
def rec_raw(H, s_, cs):
    cc = COLS[cs]; ii = H['idx'][np.ix_(s_, cc)].ravel(); m = c ** 3 + 1; Wv_ = H['Wv'][np.ix_(s_, cc)]
    return np.bincount(ii, Wv_.real.ravel(), m)[:-1] + 1j * np.bincount(ii, Wv_.imag.ravel(), m)[:-1], np.bincount(ii, H['Wd'][np.ix_(s_, cc)].ravel(), m)[:-1]
EPS = None
def rec(H, s_, cs):
    global EPS
    num, den = rec_raw(H, s_, cs)
    if EPS is None:
        d0 = rec_raw(H, np.arange(H['n']), cs)[1] * (neq0 / H['n']); EPS = np.array([np.median(d0[(sh == k) & (d0 > 0)]) if ((sh == k) & (d0 > 0)).any() else 1. for k in range(ns_)]) * a.eps
    V_ = num / (den + EPS[sh] * max(len(s_), 1) / neq0); v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
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
def crop_four(v):
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; Fv = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(Fv[np.ix_(sb, sb, sb)]))) * msk3
    return np.fft.fftn(np.fft.ifftshift(vc)).ravel()
# ---------------- arm readout on crop-c maps (arm_readout.py: lp 10, standardised, closed / open masks from the deposited closed map vs the consensus)
gq = np.fft.fftfreq(c) * c / n / a.smpd; qx, qy, qz = np.meshgrid(gq, gq, gq, indexing='ij'); GL = np.exp(-0.5 * (qx ** 2 + qy ** 2 + qz ** 2) * (10. / 1.2) ** 2).ravel()
def real_lp(Fm): return np.fft.fftshift(np.real(np.fft.ifftn((Fm * GL).reshape(c, c, c))))
cons = real_lp(crop_four(rd(a.cons))); supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def std(v): v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
consS = std(cons); clS = std(real_lp(crop_four(rd(a.closed)))); dd = clS - consS; thr = 0.5 * np.percentile(np.abs(dd[supp | (clS > 0.2 * clS.max())]), 99)
mC = dd > thr; mO = dd < -thr; refC = clS[mC].mean(); refO = consS[mO].mean()
def readout(Fm): v = std(real_lp(dock(Fm, CONS0))); return v[mC].mean() / refC, v[mO].mean() / refO
log('crop-%d arm readout: closed-arm mask %d voxels, open-arm mask %d voxels' % (c, mC.sum(), mO.sum()))
CONS0 = crop_four(rd(a.cons))
# truth occupancy per particle (synthetic): readout of each true map
occ_true = np.full(N, np.nan)
if a.truthmaps:
    TL = [l.split() for l in open(a.truthmaps) if l.strip()]; RT = np.array([readout(crop_four(rd(l[1])))[0] for l in TL]); occ_true[tt >= 0] = RT[tt[tt >= 0]]
    log('truth maps closed-site occupancy: ' + ' '.join('%d:%.2f' % (k, RT[k]) for k in range(0, len(RT), 5)))
# ---------------- held-out windows
LS = {'A': (lamA, IA), 'B': (lamB, IB)}; WIN = []
for fs, ft in (('A', 'B'), ('B', 'A')):
    ls_, Is_ = LS[fs]; lt_ = LS[ft][0]; qs = np.percentile(ls_, np.linspace(0, 100, a.nw + 1)); wl = np.clip(np.searchsorted(qs, ls_, side='right') - 1, 0, a.nw - 1)
    for w in range(a.nw):
        ro = [readout(rec(H, np.where(wl[H['pi']] == w)[0], ft)) for H in HS]; m_ = wl == w
        WIN.append(dict(f=fs, w=w, n=int(m_.sum()), ls=float(ls_[m_].mean()), lt=float(lt_[m_].mean()), occ=[r[0] for r in ro], opn=[r[1] for r in ro],
                        true=float(np.nanmean(occ_true[m_])) if np.isfinite(occ_true[m_]).any() else np.nan))
        W_ = WIN[-1]; log('  fold %s window %2d: %6d ptcls, lam_S %+.3f lam_test %+.3f | closed-site %.3f / %.3f (halves) open-site %.3f / %.3f | true closed-site %.3f' % (fs, w + 1, W_['n'],
            W_['ls'], W_['lt'], W_['occ'][0], W_['occ'][1], W_['opn'][0], W_['opn'][1], W_['true']))
# ---------------- calibration, reliability, noise scale from window means (lam oriented so that it increases with the occupancy)
sg = np.sign(np.corrcoef([np.mean(w['occ']) for w in WIN], [w['lt'] for w in WIN])[0, 1])
if sg < 0:
    lamA, lamB = -lamA, -lamB; LS = {'A': (lamA, IA), 'B': (lamB, IB)}
    for w in WIN: w['ls'], w['lt'] = -w['ls'], -w['lt']
res = dict(windows=WIN, orient=float(sg)); cr = []; pts = []
for fs in 'AB':
    Wf = [w for w in WIN if w['f'] == fs]; mo = np.array([np.mean(w['occ']) for w in Wf]); mls = np.array([w['ls'] for w in Wf]); mlt = np.array([w['lt'] for w in Wf])
    b1 = np.polyfit(mls, mo, 1)[0]; b2 = np.polyfit(mo, mlt, 1)[0]; ls_, Is_ = LS[fs]; rel = b1 * b2; c_ = ls_.var() * (1. - rel) / np.mean(1. / np.maximum(Is_, 1e-30))
    cr.append(c_); pts += list(zip(mo, mlt)); log('fold %s: b1 %.3f, b2 %.3f, reliability %.3f, c_reg %.3f' % (fs, b1, b2, rel, c_)); res['fold' + fs] = dict(b1=b1, b2=b2, rel=rel, c=c_)
c_reg = float(np.mean(cr)); pts = np.array(sorted(pts)); xs, ys = pts[:, 0], np.maximum.accumulate(pts[:, 1]); slope = np.polyfit(xs, ys, 1)[0]
g0, g1, gn = [float(x) for x in a.grid.split(':')]; og = np.linspace(g0, g1, int(gn)); lam_o = np.interp(og, xs, ys); lo, hi = og < xs[0], og > xs[-1]
lam_o[lo] = ys[0] + slope * (og[lo] - xs[0]); lam_o[hi] = ys[-1] + slope * (og[hi] - xs[-1])
lam_full = 0.5 * (lamA + lamB); I_full = IA + IB; lam_f = (lamA * IA + lamB * IB) / np.maximum(I_full, 1e-30); zsp = (lamA - lamB) / np.sqrt(1. / np.maximum(IA, 1e-30) + 1. / np.maximum(IB, 1e-30))
log('noise scale c_reg %.3f (split %.3f); calibration lam(occ) slope %.3f per unit occupancy' % (c_reg, zsp.var(), slope))
def lik(cv): dl = lam_f[:, None] - lam_o[None]; L = np.exp(-0.5 * dl ** 2 * (I_full / cv)[:, None]); return L / np.maximum(L.max(1, keepdims=True), 1e-300)
pt = None
if np.isfinite(occ_true).any():
    m_ = np.isfinite(occ_true); idx = np.clip(np.rint((occ_true[m_] - og[0]) / (og[1] - og[0])).astype(int), 0, len(og) - 1); pt = np.bincount(idx, minlength=len(og)) / m_.sum()
for vn, cv in (('c_reg', c_reg), ('c_split', float(zsp.var()))):
    fr = LD.fit_cv_boot(lik(cv), half, nboot=a.nboot, seed=a.seed); p = fr['p']; B = fr['boot']
    def mass(q, lo_, hi_): return float(q[(og >= lo_) & (og < hi_)].sum())
    row = dict(c=cv, lam=fr['lam'], p=p.tolist(), lo=np.percentile(B, 10, 0).tolist(), hi=np.percentile(B, 90, 0).tolist(),
               closed=mass(p, 0.8, 9), inter=mass(p, 0.45, 0.8), open=mass(p, -9, 0.45), closed_band=[float(np.percentile([mass(b, 0.8, 9) for b in B], q)) for q in (10, 90)],
               inter_band=[float(np.percentile([mass(b, 0.45, 0.8) for b in B], q)) for q in (10, 90)])
    res[vn] = row
    log('RESULT %s: p(occ) closed (>=0.8) %.3f [%.3f %.3f] intermediate %.3f [%.3f %.3f] open %.3f (lam_pen %.0e)%s' % (vn, row['closed'], *row['closed_band'], row['inter'], *row['inter_band'], row['open'], fr['lam'],
        '' if pt is None else ' | truth closed %.3f intermediate %.3f open %.3f, tv %.3f' % (mass(pt, 0.8, 9), mass(pt, 0.45, 0.8), mass(pt, -9, 0.45), 0.5 * np.abs(p - pt).sum())))
    log('   p(occ) on %.2f..%.2f: %s' % (og[0], og[-1], ' '.join('%.3f' % x for x in p)))
    if pt is not None: log('   truth          : %s' % ' '.join('%.3f' % x for x in pt))
res['grid'] = og.tolist(); res['lam_occ'] = lam_o.tolist(); res['c_reg'] = c_reg; res['truth'] = None if pt is None else pt.tolist()
# map-level occupancy (no deconvolution): window counts x their mean occupancy
res['map_mean_occ'] = float(np.sum([w['n'] * np.mean(w['occ']) for w in WIN]) / np.sum([w['n'] for w in WIN]))
json.dump(res, open(a.out + '/occ.json', 'w'), indent=1); log('map-level mean closed-site occupancy %.3f; DONE %s' % (res['map_mean_occ'], a.out))

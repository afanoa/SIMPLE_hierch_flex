#!/usr/bin/env python3
"""lsamp_repro_bins3.py = lsamp_repro_bins2.py for runs whose particle sets differ (subsets by state): the orientation from the
particles shared with the reference run. lsamp_repro_bins2.py = lsamp_repro_bins.py with the orientation taken from a reference run (--ref_rundir: each run oriented to correlate
positively with the reference run's full coordinate, the reference itself by its arm readout; needed for axes that are not the closed-open
direction) and a choice of the sampled axis (--axis). lsamp_repro_bins.py: the landscapes of independent tree runs of one node compared at EQUAL particle counts. Maps mode (one run per
call): the run's fold coordinates from lsamp2h.py (coords.npz: tA from pixel set A, tB from set B; oriented closed end first by the run's
final_all arm readout, as lsamp_repro_cmp.py) cut at the same quantile edges for every run; every bin reconstructed from the TEST pixels
(fold A bins from set B, fold B bins from set A: no selection imprint) -> OUT/NAME/bins_<fold>/state_NN.mrc + labels.txt, and the node
consensus of each pixel set (cons_<set>.mrc). The particle loader and reconstruction are lsamp2h.py's.
--cmp 1 (after arm_readout.py wrote OUT/NAME/arm_bins_<fold>.txt): arm readouts per bin, and whole-map agreement at equal N: deviation
maps (bin minus the consensus of the same pixels; lp 15 A, sphere) correlated across runs (run X fold A vs run Y fold B: independent
noise) against the same run's two folds; the end-to-end difference maps likewise; map-space ruler positions (direction = all runs'
closed-minus-open end bins of one fold; the other fold's maps projected; 0 / 1 = the runs' mean open / closed end bins).
usage: lsamp_repro_bins3.py --name N --rundir D --scores SC --poses P --sim .. --mr .. --cache CP --out OUT [--ref_rundir D0 --axis 1] | --cmp 1 --names A,B,.. --out OUT"""
import os, glob, argparse, itertools, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--name', default=''); ap.add_argument('--rundir', default=''); ap.add_argument('--scores', default=''); ap.add_argument('--poses', default='')
ap.add_argument('--sim', default=''); ap.add_argument('--mr', default=''); ap.add_argument('--cache', default=''); ap.add_argument('--out', required=True)
ap.add_argument('--edges', default='0,2,5,10,20,35,50,65,80,90,95,98,100', help='quantile edges (%%) from the closed end, the same for every run')
ap.add_argument('--cmp', type=int, default=0); ap.add_argument('--names', default=''); ap.add_argument('--crop', type=int, default=48)
ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.); ap.add_argument('--cs', type=float, default=2.7)
ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05); ap.add_argument('--nview', type=int, default=40)
ap.add_argument('--fsc_res', type=float, default=20.); ap.add_argument('--ref_rundir', default=''); ap.add_argument('--axis', type=int, default=1, help='sampled-axis column of coords.npz (1 = the first)')
a = ap.parse_args(); q = np.array([float(v) for v in a.edges.split(',')]); nb = len(q) - 1
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def arm(fn):   # arm_readout.py output: state number -> (closed-arm, open-arm)
    out = {}
    for l in open(fn):
        if l.startswith('#') or not l.strip(): continue
        w = l.split(); out[int(os.path.basename(w[0])[6:8])] = (float(w[1]), float(w[2]))
    return out
if a.cmp:
    names = [s for s in a.names.split(',') if s]; n = 160; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=a.smpd), indexing='ij')
    G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (15. / 1.2) ** 2); g = np.arange(n) - n // 2; sph = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2) < 0.48 * n
    def rd(p): h = np.fromfile(p, np.int32, 256); m = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=m ** 3).reshape(m, m, m).astype(np.float64)
    def lp(fn): return np.fft.irfftn(np.fft.rfftn(rd(fn)) * G, s=(n, n, n))[sph]
    def cc(x, y): x = x - x.mean(); y = y - y.mean(); return float(x @ y / np.sqrt((x @ x) * (y @ y)))
    V, AR = {}, {}
    for nm in names:
        for f_, pix in (('A', 'B'), ('B', 'A')):   # fold A bins are maps of pixel set B: their deviation is taken from the pixel-B consensus
            D = os.path.join(a.out, nm); V[nm, f_] = [lp(p) for p in sorted(glob.glob(D + '/bins_%s/state_*.mrc' % f_))]; V[nm, f_, 'c'] = lp(D + '/cons_%s.mrc' % pix)
            AR[nm, f_] = arm(D + '/arm_bins_%s.txt' % f_)
    lab0 = np.loadtxt(os.path.join(a.out, names[0], 'bins_A', 'labels.txt'), dtype=int); sz = np.bincount(lab0, minlength=nb + 1)[1:]
    log('== equal-N bins of each run\'s fold coordinates (closed end first; maps from the held-out pixels): %d bins, quantile edges %s %%, sizes %s' % (nb, a.edges, ' '.join(map(str, sz))))
    for j, w in enumerate(('closed', 'open')):
        log('-- held-out %s-arm per bin (fold A / fold B per run; across-run sd of the fold means vs the within-run fold difference)' % w)
        log('   bin  quantile   ' + '  '.join('%11s' % nm for nm in names) + '   sd(runs)  |A-B|')
        for k in range(1, nb + 1):
            m_ = [0.5 * (AR[nm, 'A'][k][j] + AR[nm, 'B'][k][j]) for nm in names]; dAB = np.mean([abs(AR[nm, 'A'][k][j] - AR[nm, 'B'][k][j]) for nm in names])
            log('   %3d  %5.1f-%-5.1f ' % (k, q[k - 1], q[k]) + '  '.join('  %.2f/%.2f' % (AR[nm, 'A'][k][j], AR[nm, 'B'][k][j]) for nm in names) + '   %.3f    %.3f' % (np.std(m_, ddof=1), dAB))
    dev = lambda nm, f_, k: V[nm, f_][k] - V[nm, f_, 'c']; prs = list(itertools.combinations(names, 2))
    log('-- deviation maps (bin minus consensus, lp 15 A): corr within a run (fold A vs fold B, mean over runs) and across runs (X fold A vs Y fold B and the reverse)')
    log('   bin  within  across(mean)  ' + '  '.join('%s/%s' % pr for pr in prs))
    for k in range(nb):
        wi = np.mean([cc(dev(nm, 'A', k), dev(nm, 'B', k)) for nm in names]); cr = [0.5 * (cc(dev(x, 'A', k), dev(y, 'B', k)) + cc(dev(x, 'B', k), dev(y, 'A', k))) for x, y in prs]
        log('   %3d  %.3f   %.3f        ' % (k + 1, wi, np.mean(cr)) + '  '.join('%.3f' % v for v in cr))
    e2e = lambda nm, f_: V[nm, f_][0] - V[nm, f_][nb - 1]
    log('-- end-to-end difference maps (closed bin minus open bin): within %s; across %s' % (' '.join('%s %.3f' % (nm, cc(e2e(nm, 'A'), e2e(nm, 'B'))) for nm in names),
        ' '.join('%s/%s %.3f' % (x, y, 0.5 * (cc(e2e(x, 'A'), e2e(y, 'B')) + cc(e2e(x, 'B'), e2e(y, 'A')))) for x, y in prs)))
    P_ = {}
    for fd, fp in (('B', 'A'), ('A', 'B')):   # ruler from fold fd's end bins (all runs), positions of fold fp's maps (the other pixel set)
        d = sum(e2e(nm, fd) for nm in names); pr_ = {nm: np.array([V[nm, fp][k] @ d for k in range(nb)]) for nm in names}
        lo, hi = np.mean([pr_[nm][nb - 1] for nm in names]), np.mean([pr_[nm][0] for nm in names])
        for nm in names: P_[nm, fp] = (pr_[nm] - lo) / (hi - lo)
    log('-- map-space ruler (whole-map closed-open direction; 0 = mean open end bin, 1 = mean closed end bin; mean of the two fold directions)')
    log('   bin  ' + '  '.join('%7s' % nm for nm in names) + '   sd(runs)  |foldA-foldB|')
    for k in range(nb):
        pk = [0.5 * (P_[nm, 'A'][k] + P_[nm, 'B'][k]) for nm in names]
        log('   %3d  ' % (k + 1) + '  '.join('%7.3f' % v for v in pk) + '   %.3f     %.3f' % (np.std(pk, ddof=1), np.mean([abs(P_[nm, 'A'][k] - P_[nm, 'B'][k]) for nm in names])))
    raise SystemExit
# ---------------- one run: node rows and fold coordinates (coords.npz), oriented closed end first
M = np.load(a.scores + '/meta.npz'); NM = len(M['eo']); C = np.load(a.rundir + '/coords.npz'); rows = C['rows']; N = len(rows); eo = M['eo'][rows]
def orient(d):   # +1 / -1: the run's first sample is its closed end (final_all arm readout)
    af = arm(d + '/arm_final_all.txt'); return -1. if af[1][0] < af[max(af)][0] else 1.
j = a.axis - 1
if a.ref_rundir and os.path.abspath(a.ref_rundir) != os.path.abspath(a.rundir):
    C0 = np.load(a.ref_rundir + '/coords.npz'); sh_, i1, i0 = np.intersect1d(rows, C0['rows'], return_indices=True)   # shared particles orient the run
    sg = 1. if np.corrcoef(C['tfull'][i1, j], orient(a.ref_rundir) * C0['tfull'][i0, min(j, C0['tfull'].shape[1] - 1)])[0, 1] >= 0 else -1.
else: sg = orient(a.rundir)
T = {f_: sg * C['t%s' % f_][:, j] for f_ in 'AB'}; log('%s: %d particles, %d sampled axes (column %d used), orientation %+d' % (a.name, N, C['tA'].shape[1], a.axis, int(sg)))
# ---------------- particle data (crop 48, the scoring cache): halves = half stacks or the eo partition
lamb = 12.2643247 / np.sqrt(a.kv * 1e3 * (1. + a.kv * 1e3 * 0.978466e-6))
def rotm(e1, e2, e3):
    a1, a2, a3 = np.radians(e1), np.radians(e2), np.radians(e3); c1, si1, c2, si2, c3, si3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Z1 = np.zeros((len(e1), 3, 3)); Z1[:, 0, 0] = c1; Z1[:, 0, 1] = -si1; Z1[:, 1, 0] = si1; Z1[:, 1, 1] = c1; Z1[:, 2, 2] = 1
    Y2 = np.zeros((len(e1), 3, 3)); Y2[:, 0, 0] = c2; Y2[:, 0, 2] = si2; Y2[:, 2, 0] = -si2; Y2[:, 2, 2] = c2; Y2[:, 1, 1] = 1
    Z3 = np.zeros((len(e1), 3, 3)); Z3[:, 0, 0] = c3; Z3[:, 0, 1] = -si3; Z3[:, 1, 0] = si3; Z3[:, 1, 1] = c3; Z3[:, 2, 2] = 1
    return Z1 @ Y2 @ Z3
SIMS = a.sim.split(','); full = a.mr == 'full'; n = int(np.frombuffer(open(SIMS[0] + '/all.mrcs', 'rb').read(1024), np.int32)[0]); c = a.crop
R0 = a.mskdiam / 2 / a.smpd; fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); kk0 = np.stack([fx2, fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
r2d = np.hypot(kk0[:, 0], kk0[:, 1]); pcols = np.where(r2d <= c // 2 - 0.5)[0]; pr = r2d[pcols]; psh = np.round(pr).astype(int)
ppar = np.mod(np.round(kk0[pcols, 0]).astype(int) + np.round(kk0[pcols, 1]).astype(int), 2) == 1   # True = pixel set B
f3 = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(f3, f3, f3, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx ** 2 + ky ** 2 + kz ** 2) * n).astype(int), c // 2).ravel()
g3 = np.arange(c) - c // 2; r3 = np.sqrt(g3[:, None, None] ** 2 + g3[None, :, None] ** 2 + g3[None, None, :] ** 2) * n / c; msk3 = np.clip((R0 + 4 - r3) / 8., 0, 1).astype(np.float32)
ns_ = c // 2 + 1; khi = c // 2 - 1; kres = int(round(n * a.smpd / a.fsc_res)); KV = np.stack([np.arange(c ** 3) % c, (np.arange(c ** 3) // c) % c, np.arange(c ** 3) // (c * c)], 1); KV = np.where(KV >= c // 2, KV - c, KV).astype(np.float64) / n
FB = (sh >= 2) & (sh <= khi)
P_all = np.array([[float(kv(l).get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')] for l in open(a.poses)]); assert len(P_all) == NM
inrow = np.full(NM, -1); inrow[rows] = np.arange(N)
def load_half(tag, crow, mrow):
    Fc = np.load('%s_%s_c%d.npy' % (a.cache, tag, c), mmap_mode='r'); pp = P_all[mrow]; nn = len(crow)
    Wv = np.zeros((nn, len(pcols)), np.complex64); Wd = np.zeros((nn, len(pcols)), np.float32); idx = np.zeros((nn, len(pcols)), np.int32); pw = np.zeros((2, ns_)); pc = np.zeros((2, ns_))
    kf = (np.hypot(fx2, fy2) / a.smpd).ravel()[pcols]; ang = np.arctan2(fy2, fx2).ravel()[pcols]
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
    ph, th = np.radians(pp[:, 0]), np.radians(pp[:, 1]); v = np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
    i = np.arange(a.nview) + 0.5; z = 1 - 2 * i / a.nview; rr = np.sqrt(np.maximum(1 - z * z, 0)); gg = np.pi * (3 - np.sqrt(5)) * i
    vb = np.argmax(v @ np.c_[rr * np.cos(gg), rr * np.sin(gg), z].T, 1); dfm = 0.5 * (pp[:, 5] + pp[:, 6]); dq = np.searchsorted(np.quantile(dfm, [1 / 3., 2 / 3.]), dfm)
    return dict(Wv=Wv, Wd=Wd, idx=idx, ti=inrow[mrow], st=vb * 3 + dq, n=nn, w=(1. / np.maximum(pw / np.maximum(pc, 1), 1e-30)).astype(np.float32))
if full: HS = [load_half('A', rows[eo == e], rows[eo == e]) for e in (0, 1)]
else:
    HS = []
    for h, tg in enumerate('AB'):
        mr = np.loadtxt(a.mr.split(',')[h]).astype(int) - 1; sel = np.where(inrow[mr] >= 0)[0]; HS.append(load_half(tg, sel, mr[sel]))
log('map data: halves %d / %d particles (%s), crop %d, %d pixels per image (%d set A / %d set B)' % (HS[0]['n'], HS[1]['n'], 'eo partition' if full else 'half stacks', c, len(pcols), (~ppar).sum(), ppar.sum()))
COLS = {'all': np.arange(len(pcols)), 'A': np.where(~ppar)[0], 'B': np.where(ppar)[0]}; neq0 = 5000
def rec_raw(H, sel, cs):
    ii = H['idx'][np.ix_(sel, cs)].ravel(); m = c ** 3 + 1; Wv_ = H['Wv'][np.ix_(sel, cs)]
    return np.bincount(ii, Wv_.real.ravel(), m)[:-1] + 1j * np.bincount(ii, Wv_.imag.ravel(), m)[:-1], np.bincount(ii, H['Wd'][np.ix_(sel, cs)].ravel(), m)[:-1]
def shell_eps(H):
    den = rec_raw(H, np.arange(H['n']), COLS['all'])[1] * (neq0 / H['n']); return np.array([np.median(den[(sh == k) & (den > 0)]) if ((sh == k) & (den > 0)).any() else 1. for k in range(ns_)]) * a.eps
EPS = 0.5 * (shell_eps(HS[0]) + shell_eps(HS[1]))
def rec(H, sel, cs='all'):   # map of particles sel (indices into the half) from the pixels of set cs, fixed regulariser x N / N0, masked in real space
    cs_ = COLS[cs]; num, den = rec_raw(H, sel, cs_); f_ = len(cs_) / len(pcols); V_ = num / (den + EPS[sh] * f_ * max(len(sel), 1) / neq0)
    v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel()
CONS = {cs: [rec(H, np.arange(H['n']), cs) for H in HS] for cs in ('all', 'A', 'B')}
def write_mrc(fn, v, apix):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); h[0:3] = [nx, ny, nz]; h[3] = 2; h[7:10] = [nx, ny, nz]
    hf = h.view(np.float32); hf[10:13] = [nx * apix, ny * apix, nz * apix]; hf[13:16] = 90.; h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = v.min(), v.max(), v.mean()
    h[52] = np.frombuffer(b'MAP ', np.int32)[0]; h[53] = np.frombuffer(bytes([0x44, 0x44, 0, 0]), np.int32)[0]
    with open(fn, 'wb') as f: f.write(h.tobytes()); f.write(v.tobytes())
NB = 160
def padmap(M_):
    Fc = M_.reshape(c, c, c); Fb = np.zeros((NB, NB, NB), np.complex128); i_ = np.r_[0:c // 2, NB - c // 2:NB]; Fb[np.ix_(i_, i_, i_)] = Fc
    return np.fft.fftshift(np.real(np.fft.ifftn(Fb))) * (NB / c) ** 3
# ---------------- fixed-quantile bins of each fold coordinate, maps from the test pixels
D0 = os.path.join(a.out, a.name); os.makedirs(D0, exist_ok=True)
for cs_ in ('A', 'B'): write_mrc(os.path.join(D0, 'cons_%s.mrc' % cs_), padmap(0.5 * (CONS[cs_][0] + CONS[cs_][1])), n * a.smpd / NB)
for f_, test in (('A', 'B'), ('B', 'A')):
    t = T[f_]; lab = np.searchsorted(np.percentile(t, q[1:-1]), t, side='right') + 1; D = os.path.join(D0, 'bins_%s' % f_); os.makedirs(D, exist_ok=True)
    lo_ = np.zeros(NM, int); lo_[rows] = lab; np.savetxt(os.path.join(D, 'labels.txt'), lo_, fmt='%d'); labh = [lab[H['ti']] for H in HS]
    for k in range(1, nb + 1):
        M_ = 0.5 * (rec(HS[0], np.where(labh[0] == k)[0], test) + rec(HS[1], np.where(labh[1] == k)[0], test)); write_mrc(os.path.join(D, 'state_%02d.mrc' % k), padmap(M_), n * a.smpd / NB)
    log('%s fold %s: %d bins of coordinate %s, maps from pixels %s, sizes %s' % (a.name, f_, nb, f_, test, ' '.join(str((lab == k).sum()) for k in range(1, nb + 1))))
log('DONE')

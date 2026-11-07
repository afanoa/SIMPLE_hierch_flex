#!/usr/bin/env python3
"""lsamp_polemaps.py: pole maps for a run's landscape score space. Per pole (a soft site ranking's posterior per master row, e.g. the
factorial session's <run>_z1lo/post_master.txt), the top-N particles of the node, reconstructed at the RUN's poses (the score space's pose
table) from the scoring cache (crop 48 = the score band), even map from eo 0 and odd map from eo 1 of the run's eo flags (the score's
cross-half convention), each pole's per-shell signal amplitude (even x odd cross power) matched to the mean of the run's leaf maps (one
global per-shell scale serves every reference in lsamp_score2.py), written at the stack box. OUT gets the pole maps, leaves_list.txt = the
run's leaves + the poles (PREFIX + name, so lsamp2h.py --use PREFIX takes them) and final_labels.txt = the run's labels (poles own no rows).
The particle loader and reconstruction are lsamp2h.py's.
usage: lsamp_polemaps.py --leaves LIST --labels LAB --poses P --sim .. --mr .. --cache CP --poles NAME:POST:N,.. --prefix r_ --out OUT"""
import os, argparse, numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--leaves', required=True); ap.add_argument('--labels', required=True); ap.add_argument('--poses', required=True); ap.add_argument('--poles', required=True)
ap.add_argument('--prefix', required=True); ap.add_argument('--out', required=True); ap.add_argument('--sim', required=True); ap.add_argument('--mr', required=True); ap.add_argument('--cache', required=True)
ap.add_argument('--crop', type=int, default=48); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--kv', type=float, default=300.)
ap.add_argument('--cs', type=float, default=2.7); ap.add_argument('--fraca', type=float, default=0.1); ap.add_argument('--ctfabs', type=int, default=0); ap.add_argument('--eps', type=float, default=0.05)
ap.add_argument('--nview', type=int, default=40); ap.add_argument('--fsc_res', type=float, default=20.)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def rd(p): h = np.fromfile(p, np.int32, 256); m = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=m ** 3).reshape(m, m, m).astype(np.float64)
LAB = np.loadtxt(a.labels, dtype=int); NM = len(LAB); rows = np.where(LAB > 0)[0]; N = len(rows)
EOm = np.array([int(float(kv(l).get('eo', -1))) for l in open(a.poses)]); eo = EOm[rows]
log('node: %d particles (eo 0 / 1: %d / %d)' % (N, (eo == 0).sum(), (eo == 1).sum()))
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
def write_mrc(fn, v, apix):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); h[0:3] = [nx, ny, nz]; h[3] = 2; h[7:10] = [nx, ny, nz]
    hf = h.view(np.float32); hf[10:13] = [nx * apix, ny * apix, nz * apix]; hf[13:16] = 90.; h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = v.min(), v.max(), v.mean()
    h[52] = np.frombuffer(b'MAP ', np.int32)[0]; h[53] = np.frombuffer(bytes([0x44, 0x44, 0, 0]), np.int32)[0]
    with open(fn, 'wb') as f: f.write(h.tobytes()); f.write(v.tobytes())
NB = n
def padmap(M_):
    Fc = M_.reshape(c, c, c); Fb = np.zeros((NB, NB, NB), np.complex128); i_ = np.r_[0:c // 2, NB - c // 2:NB]; Fb[np.ix_(i_, i_, i_)] = Fc
    return np.fft.fftshift(np.real(np.fft.ifftn(Fb))) * (NB / c) ** 3
# ---------------- pole maps: top-N of each posterior, per eo half of the run's flags, amplitude-matched to the leaves
def rec_eo(e, flag):   # flag over the node rows; both half stacks pooled (a run's eo split need not follow the stacks)
    num = 0.; den = 0.; nsel = 0
    for H in HS:
        sel = np.where(flag[H['ti']] & (eo[H['ti']] == e))[0]; nu, de = rec_raw(H, sel, COLS['all']); num = num + nu; den = den + de; nsel += len(sel)
    V_ = num / (den + EPS[sh] * nsel / neq0); v = np.fft.fftshift(np.real(np.fft.ifftn(V_.reshape(c, c, c)))) * msk3; return np.fft.fftn(np.fft.ifftshift(v)).ravel(), nsel
def crop4(v):   # lsamp_score2.py's crop_four without its padding element
    nb = v.shape[0]; sb = np.r_[0:c // 2, nb - c // 2:nb]; F = np.fft.fftn(np.fft.ifftshift(v)); vc = np.fft.fftshift(np.real(np.fft.ifftn(F[np.ix_(sb, sb, sb)]))); return np.fft.fftn(np.fft.ifftshift(vc)).ravel()
def xpow(E, O): return np.bincount(sh, np.real(E * np.conj(O)), ns_)
L = [l.split() for l in open(a.leaves) if l.strip() and not l.startswith('#')]; nb = rd(L[0][1]).shape[0]
SL = np.mean([xpow(crop4(rd(l[1])), crop4(rd(l[2]))) for l in L], 0); kb = int(round(n * a.smpd / 12.))
log('%d leaf maps (box %d): signal power per shell from their even x odd cross power; band to 12 A = shell %d' % (len(L), nb, kb))
out = [' '.join(l) for l in L]
for spec in a.poles.split(','):
    nm, post, ntop = spec.split(':'); p = np.loadtxt(post)[rows]; ntop = int(ntop); flag = np.zeros(N, bool); flag[np.argsort(-p)[:ntop]] = True
    (E, ne), (O, no) = rec_eo(0, flag), rec_eo(1, flag); SP = xpow(E, O); ok = (SL > 0) & (SP > 0); f = np.where(ok, np.sqrt(np.where(ok, SL, 1.) / np.where(ok, SP, 1.)), 1.)
    log('pole %s: top %d of %s (eo 0 / 1: %d / %d); amplitude match per shell 2..%d: %s' % (nm, ntop, post, ne, no, kb, ' '.join('%.2f' % v for v in f[2:kb + 1])))
    fn = {}
    for tg, M_ in (('even', E), ('odd', O)):   # crop4(file) returns the matched crop map: padmap scales by (NB/c)^3, undone here
        fn[tg] = os.path.join(a.out, '%s%s_%s.mrc' % (a.prefix, nm, tg)); write_mrc(fn[tg], padmap(M_ * f[sh]) / (NB / c) ** 3, n * a.smpd / NB)
    out.append('%s%s %s %s %d' % (a.prefix, nm, fn['even'], fn['odd'], ntop))
open(os.path.join(a.out, 'leaves_list.txt'), 'w').write('\n'.join(out) + '\n'); np.savetxt(os.path.join(a.out, 'final_labels.txt'), LAB, fmt='%d')
log('wrote %s: %d leaves + %d poles; DONE' % (os.path.join(a.out, 'leaves_list.txt'), len(L), len(out) - len(L)))

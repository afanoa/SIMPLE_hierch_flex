#!/usr/bin/env python3
"""land_morph.py: a PHYSICAL path between two maps by unbalanced optimal-transport displacement interpolation (no bodies or masks drawn by
hand). Both maps are low-passed to LP, made non-negative and B is scaled to A on their shared support; the shared density is min(A, B);
the moving parts dA = A - shared and dB = B - shared are point clouds (voxels of a STEP-px grid above FRAC of their maximum) matched by
entropic unbalanced OT (log-domain Sinkhorn, blur EPS A^2, marginal relaxation REACH A: mass with no partner within ~REACH fades
instead of travelling). Step s places each transported element at (1-s) x_i + s y_j (McCann interpolation), splats it trilinearly and
blurs by STEP/2 px; untransported mass of A fades as (1-s), of B appears as s. path(s) = shared + moving(s). Also writes the linear
cross-fade (1-s) A + s B (what a two-state mixture looks like) for comparison.
With --truth DIR (conf<ID0+k>_*.mrc, k = 0..NSTEPS-1, the same path sampled at equal steps): per step the correlation of the OT path and
of the cross-fade with the truth inside the moving region, at lp LP.
usage: land_morph.py A B OUT [--nsteps 31 --lp 8 --step 2 --frac 0.05 --eps 2 --reach 40 --truth DIR --id0 200 --mask MSK]"""
import os, sys, glob, argparse, numpy as np
from scipy.ndimage import gaussian_filter, binary_dilation
from scipy.special import logsumexp
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, write_mrc
ap = argparse.ArgumentParser(); ap.add_argument('A'); ap.add_argument('B'); ap.add_argument('out'); ap.add_argument('--nsteps', type=int, default=31)
ap.add_argument('--lp', type=float, default=8.); ap.add_argument('--step', type=int, default=2); ap.add_argument('--frac', type=float, default=0.05)
ap.add_argument('--eps', type=float, default=2.); ap.add_argument('--reach', type=float, default=40.); ap.add_argument('--truth', default=''); ap.add_argument('--id0', type=int, default=200)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--nit', type=int, default=400); ap.add_argument('--maxpts', type=int, default=8000)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def log(*s): print(*s, flush=True)
A0 = read_mrc(a.A); B0 = read_mrc(a.B); n = A0.shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij')
G = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (a.lp / 1.2) ** 2); lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G))
A = np.maximum(lp(A0), 0); B = np.maximum(lp(B0), 0)
sup = (A > 0.2 * A.max()) & (B > 0.2 * B.max()); B *= (A[sup] @ B[sup]) / max(B[sup] @ B[sup], 1e-30)
S = np.minimum(A, B); dA = A - S; dB = B - S
log('maps %s -> %s: masses A %.4g B %.4g, shared %.4g, moving A %.4g (%.1f %%) B %.4g (%.1f %%)' % (os.path.basename(a.A), os.path.basename(a.B), A.sum(), B.sum(), S.sum(), dA.sum(),
    100 * dA.sum() / A.sum(), dB.sum(), 100 * dB.sum() / B.sum()))
def cloud(D):   # coarse-grid point cloud (A) and weights of a non-negative density
    st = a.step; Dc = D.reshape(n // st, st, n // st, st, n // st, st).sum((1, 3, 5)); thr = a.frac * Dc.max()
    if (Dc > thr).sum() > a.maxpts: thr = np.sort(Dc.ravel())[-a.maxpts]   # cap the cloud at the MAXPTS heaviest cells
    m = Dc > thr; idx = np.argwhere(m)
    return (idx * st + (st - 1) / 2.) * a.smpd, Dc[m], m
X, wa, _ = cloud(dA); Y, wb, _ = cloud(dB); ma, mb = wa.sum(), wb.sum(); pa, pb = wa / ma, wb / mb
log('point clouds: %d (A) / %d (B) at %d-px grid, %.1f / %.1f %% of the moving mass' % (len(X), len(Y), a.step, 100 * ma / dA.sum(), 100 * mb / dB.sum()))
C = np.maximum((X ** 2).sum(1)[:, None] + (Y ** 2).sum(1)[None, :] - 2. * X @ Y.T, 0.); rho = a.reach ** 2; la, lb = np.log(pa), np.log(pb)
fa = np.zeros(len(X)); gb = np.zeros(len(Y)); e0 = max(C.max() / 10., a.eps ** 2)
for it in range(a.nit):   # unbalanced Sinkhorn (KL marginals), log domain, with blur annealed to EPS^2
    eps = max(a.eps ** 2, e0 * 0.97 ** it); kap = rho / (rho + eps)
    fa = -kap * eps * logsumexp(lb[None, :] + (gb[None, :] - C) / eps, axis=1)
    gb = -kap * eps * logsumexp(la[:, None] + (fa[:, None] - C) / eps, axis=0)
P = np.exp(la[:, None] + lb[None, :] + (fa[:, None] + gb[None, :] - C) / eps); Pa, Pb = P.sum(1), P.sum(0)
d = np.sqrt(C); tm = P.sum(); md = (P * d).sum() / tm; o = np.argsort(d.ravel()); cq = np.cumsum(P.ravel()[o]) / tm
log('OT plan (final blur %.1f A): transported mass %.3f of A / %.3f of B; mean displacement %.1f A (90 %% below %.1f A)' % (np.sqrt(eps), Pa.sum(), Pb.sum(), md,
    d.ravel()[o][np.searchsorted(cq, 0.9)]))
keep = P > 1e-4 * P.max(); ii, jj = np.nonzero(keep); w = P[ii, jj]
ra = np.clip(pa - Pa, 0, None); rb = np.clip(pb - Pb, 0, None)   # untransported mass (fades / appears)
def splat(pts, wts):
    g = pts / a.smpd; i0 = np.floor(g).astype(int); fr = g - i0; out = np.zeros((n, n, n))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                ww = wts * (fr[:, 0] if dx else 1 - fr[:, 0]) * (fr[:, 1] if dy else 1 - fr[:, 1]) * (fr[:, 2] if dz else 1 - fr[:, 2]); ix = np.clip(i0 + [dx, dy, dz], 0, n - 1)
                np.add.at(out, (ix[:, 0], ix[:, 1], ix[:, 2]), ww)
    return gaussian_filter(out, a.step / 2.)
MA, MB = dA.sum(), dB.sum(); fadeA = splat(X, ra) * MA; fadeB = splat(Y, rb) * MB   # untransported mass, as fractions of the moving masses
lines = []; tr = []
if a.truth:
    TF = [glob.glob('%s/conf%d_*.mrc' % (a.truth, a.id0 + k))[0] for k in range(a.nsteps)]; T0 = [np.maximum(lp(read_mrc(t)), 0) for t in TF]
    dev = np.abs(np.array(T0) - T0[0]).max(0); R = binary_dilation(dev > 0.1 * dev.max(), iterations=2)
    log('truth: %d maps, moving region %d voxels' % (len(TF), R.sum()))
def corr(u, v, m): u = u[m] - u[m].mean(); v = v[m] - v[m].mean(); return float(u @ v / np.sqrt((u @ u) * (v @ v) + 1e-30))
for k in range(a.nsteps):
    s = k / (a.nsteps - 1.); pos = (1 - s) * X[ii] + s * Y[jj]; mv = splat(pos, w) * ((1 - s) * MA + s * MB)
    Pm = S + mv + (1 - s) * fadeA + s * fadeB; Xf = (1 - s) * A + s * B
    write_mrc('%s/path_%02d.mrc' % (a.out, k), Pm.astype(np.float32), a.smpd); write_mrc('%s/xfade_%02d.mrc' % (a.out, k), Xf.astype(np.float32), a.smpd)
    if a.truth:
        c_ot, c_xf = corr(Pm, T0[k], R), corr(Xf, T0[k], R); tr.append((c_ot, c_xf)); lines.append('%d %.4f %.4f %.4f' % (k, s, c_ot, c_xf))
    else: lines.append('%d %.4f' % (k, s))
open(a.out + '/path.txt', 'w').write('# step s' + (' corr_OT corr_crossfade (truth, moving region, lp %.0f)' % a.lp if a.truth else '') + '\n' + '\n'.join(lines) + '\n')
if a.truth:
    tr = np.array(tr); log('truth correlation in the moving region per step (OT / cross-fade): ' + ' '.join('%d:%.3f/%.3f' % (k, tr[k, 0], tr[k, 1]) for k in range(0, a.nsteps, 3)))
    log('RESULT mean over the interior steps: OT %.3f, cross-fade %.3f; minimum OT %.3f, cross-fade %.3f' % (tr[1:-1, 0].mean(), tr[1:-1, 1].mean(), tr[1:-1, 0].min(), tr[1:-1, 1].min()))
log('wrote %d path maps (+ cross-fades) to %s; DONE' % (a.nsteps, a.out))

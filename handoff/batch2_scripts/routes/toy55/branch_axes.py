#!/usr/bin/env python3
"""Axes of one branch of states and even sampling along them.
The branch's states (even/odd half-maps, particle counts) -> the particle-weighted, centred cross-half Gram matrix of the maps inside
MASK over the band [HP, LP] (Re <E_s - Ebar, O_t - Obar> is unbiased for the true maps' inner products, no noise on the diagonal) ->
eigen-axes (classical MDS of the unbiased distances); per state its coordinate on every axis; per axis a mode map (the particle-weighted
regression of the centred full maps on the coordinate, MRC). Even sampling along AXIS: the coordinate range from the lowest to the
highest state cut into K equal-width bins, each bin's states pooled; the end bins hold the extreme states. Written per K as
OUT/axis<A>_k<K>/{groups.txt, states.txt} (per-row labels for mu_rec.sbatch).
usage: branch_axes.py OUT --rdir HALVES --states 6,24,... --groups groups.txt --labels states.txt --mask MASK
       [--arm arm_states.txt --hp 40 --lp 10 --k 4,5,6 --axis 1 --smpd 0]"""
import os, sys, re, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rdir', required=True); ap.add_argument('--states', required=True)
ap.add_argument('--groups', required=True); ap.add_argument('--labels', required=True); ap.add_argument('--mask', required=True); ap.add_argument('--arm', default='')
ap.add_argument('--hp', type=float, default=40.); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--k', default='4,5,6'); ap.add_argument('--axis', type=int, default=1)
ap.add_argument('--smpd', type=float, default=0.); ap.add_argument('--grid', type=int, default=0, help='G: even G x G grid over axes 1-2 (the extremes of both axes kept)'); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, v, smpd):
    n = v.shape[0]; h = np.zeros(256, np.int32); h[0:3] = n; h[3] = 2; h[7:10] = n; hf = h.view(np.float32); hf[10:13] = n * smpd; hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19], hf[20], hf[21] = v.min(), v.max(), v.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x00004444; open(p, 'wb').write(h.tobytes() + v.astype(np.float32).tobytes())
S = [s for s in a.states.split(',') if s]; cnt = {}
for l in open(a.groups):
    if l.strip() and not l.startswith('#'): w = l.split(); cnt[w[0]] = float(w[1])
arm = {}
for l in open(a.arm) if a.arm else []:
    if l.strip() and not l.startswith('#'): w = l.split(); arm[str(int(re.findall(r'(\d+)\.mrc$', w[0])[0]))] = float(w[1])
n_ = np.array([cnt[s] for s in S]); wt = n_ / n_.sum(); ev = ['%s/recvol_state%02d_even.mrc' % (a.rdir, int(s)) for s in S]; od = [e.replace('_even.mrc', '_odd.mrc') for e in ev]
E0 = rd(ev[0]); n = E0.shape[0]; smpd = a.smpd if a.smpd > 0 else float(np.fromfile(ev[0], np.float32, 13)[10]) / n
M = np.clip(rd(a.mask), 0, 1); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
res = np.where(kr > 0, smpd / np.maximum(kr, 1e-9), 1e9); band = (res <= a.hp) & (res >= a.lp); wgt = np.where((kz[band] == 0) | (np.abs(kz[band]) == 0.5), 1., 2.)
print('%d states (%d particles), box %d, smpd %.3f, band %.0f-%.0f A, mask volume fraction %.3f' % (len(S), n_.sum(), n, smpd, a.hp, a.lp, M.mean()))
Ev = np.array([np.fft.rfftn(M * rd(p))[band] for p in ev]); Ov = np.array([np.fft.rfftn(M * rd(p))[band] for p in od])
Ec, Oc = Ev - wt @ Ev, Ov - wt @ Ov; G = np.real((Ec * wgt) @ Oc.conj().T); G = (G + G.T) / 2   # unbiased centred Gram of the true maps
sw = np.sqrt(wt); A = sw[:, None] * G * sw[None, :]; lam, U = np.linalg.eigh(A); lam, U = lam[::-1], U[:, ::-1]
pos = lam[lam > 0].sum(); X = U * np.sqrt(np.maximum(lam, 0))[None, :] / sw[:, None]   # sum_s w_s X[s,k]^2 = lambda_k
v = np.array([arm.get(s, np.nan) for s in S]); rk = lambda x: np.argsort(np.argsort(x)).astype(float)
print('axis variance shares: %s (negative eigenvalues = noise: %s)' % (' '.join('%.2f' % (l / pos) for l in lam[:6]), ' '.join('%.2f' % (l / pos) for l in lam if l < 0)))
for k in range(min(4, len(S))):
    if arm and np.isfinite(v).all() and np.corrcoef(X[:, k], v)[0, 1] < 0: X[:, k] *= -1; U[:, k] *= -1   # display only: orient so the closed arm grows along +axis
    print('axis %d: Spearman with the closed-arm readout %+.2f' % (k + 1, np.corrcoef(rk(X[:, k]), rk(v))[0, 1] if arm else float('nan')))
ax = a.axis - 1; o = np.argsort(X[:, ax])
print('states along axis %d (state, coordinate, axis 2, particles, closed-arm):' % a.axis)
for i in o: print('   state_%02d  %+8.3f  %+8.3f  %6d  %.2f' % (int(S[i]), X[i, ax] / np.sqrt(max(lam[ax], 1e-30)), X[i, 1] / np.sqrt(max(lam[1], 1e-30)) if len(S) > 1 else 0, n_[i], v[i]))
np.savetxt(a.out + '/coords.txt', np.c_[[int(s) for s in S], n_, X[:, :min(4, len(S))], v], fmt='%g', header='state particles axis1..axis4 (sum_s w_s x^2 = eigenvalue) closed-arm')
# mode maps: the regression of the centred full maps on each coordinate (real space, whole box), in units of one sd of the coordinate
F = np.array([(rd(e) + rd(o_)) / 2 for e, o_ in zip(ev, od)]); Fc = F - np.tensordot(wt, F, 1)
for k in range(min(3, len(S))):
    sd = np.sqrt(max(lam[k], 1e-30)); mode = np.tensordot(wt * X[:, k] / sd, Fc, 1); wr('%s/axis%d_mode.mrc' % (a.out, k + 1), mode, smpd)
    wr('%s/axis%d_minus.mrc' % (a.out, k + 1), np.tensordot(wt, F, 1) - 2 * mode, smpd); wr('%s/axis%d_plus.mrc' % (a.out, k + 1), np.tensordot(wt, F, 1) + 2 * mode, smpd)
# even sampling along the chosen axis
lab = np.loadtxt(a.labels, dtype=int); x = X[:, ax]
for K in [int(k) for k in a.k.split(',') if k]:
    edges = np.linspace(x.min(), x.max(), K + 1); b = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, K - 1)
    grp = [[S[i] for i in np.where(b == j)[0]] for j in range(K)]; grp = [g for g in grp if g]
    d = '%s/axis%d_k%d' % (a.out, a.axis, K); os.makedirs(d, exist_ok=True); lut = np.zeros(lab.max() + 1, int)
    print('== axis %d in %d equal-width bins -> %d samples (low -> high):' % (a.axis, K, len(grp)))
    with open(d + '/groups.txt', 'w') as fo:
        fo.write('# even sample along axis %d of the branch (K %d): particles, members (merge-up states)\n' % (a.axis, K))
        for gi, g in enumerate(grp, 1):
            ix = [S.index(s) for s in g]; ng = n_[ix].sum(); fo.write('%d %d %s\n' % (gi, ng, '+'.join(g))); lut[[int(s) for s in g]] = gi
            print('   sample %d: states %s, %d particles, coordinate %.2f..%.2f, closed-arm (particle-weighted states) %.2f' % (gi, ' '.join(g), ng, x[ix].min() / np.sqrt(lam[ax]), x[ix].max() / np.sqrt(lam[ax]), (n_[ix] * v[ix]).sum() / ng))
    np.savetxt(d + '/states.txt', lut[np.clip(lab, 0, None)] * (lab > 0), fmt='%d')
# even sampling on a G x G grid over axes 1 and 2: equal-width cells between the extreme states of each axis, non-empty cells pooled
if a.grid > 0:
    G_ = a.grid; z1, z2 = X[:, 0] / np.sqrt(lam[0]), X[:, 1] / np.sqrt(lam[1]); c1 = np.clip(((z1 - z1.min()) / (np.ptp(z1) + 1e-12) * G_).astype(int), 0, G_ - 1)
    c2 = np.clip(((z2 - z2.min()) / (np.ptp(z2) + 1e-12) * G_).astype(int), 0, G_ - 1); cells = sorted(set(zip(c1, c2)))
    d = '%s/grid12_g%d' % (a.out, G_); os.makedirs(d, exist_ok=True); lut = np.zeros(lab.max() + 1, int)
    print('== %d x %d grid over axes 1-2 -> %d samples (cell axis1, axis2):' % (G_, G_, len(cells)))
    with open(d + '/groups.txt', 'w') as fo:
        fo.write('# even %d x %d grid over axes 1-2 of the branch: particles, members (merge-up states)\n' % (G_, G_))
        for gi, cc in enumerate(cells, 1):
            ix = [i for i in range(len(S)) if (c1[i], c2[i]) == cc]; g = [S[i] for i in ix]; ng = n_[ix].sum(); fo.write('%d %d %s\n' % (gi, ng, '+'.join(g))); lut[[int(x_) for x_ in g]] = gi
            print('   sample %d: cell (%d, %d), states %s, %d particles, closed-arm (particle-weighted states) %.2f' % (gi, cc[0] + 1, cc[1] + 1, ' '.join(g), ng, (n_[ix] * v[ix]).sum() / ng))
    np.savetxt(d + '/states.txt', lut[np.clip(lab, 0, None)] * (lab > 0), fmt='%d')

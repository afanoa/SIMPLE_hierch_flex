#!/usr/bin/env python3
"""3D states from direction-bin sub-class averages with KNOWN projection directions (no orientation search): a K-state
EM over the averages. M-step = Fourier-slice insertion of each state's member averages (member-weighted) into a coarse
3D grid, one grid per half (even / odd averages); E-step = cross-half: every sub-class is scored against the slices of
the other half's state volumes at its own direction (band-limited correlation inside a soft mask) and assigned to the
best state. Particles inherit their sub-class's state. The Euler/slice convention is chosen by matching the consensus
projections to the full averages.
usage: avgs_states3d.py STITCH_DIR CLSEXP_DIR CONSENSUS OUT --k 6 [--crop 80 --lp 12 --hp 40 --iters 20 --seed 0]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
ap = argparse.ArgumentParser(); ap.add_argument('stitch'); ap.add_argument('clsexp'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--k', type=int, default=6)
ap.add_argument('--crop', type=int, default=80); ap.add_argument('--lp', type=float, default=12.); ap.add_argument('--hp', type=float, default=60.); ap.add_argument('--iters', type=int, default=20); ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--npart', type=int, default=0, help='pad the particle label file to this many rows'); a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
def read_stack(fn):
    h = np.frombuffer(open(fn, 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
    return np.array(np.memmap(fn, np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)), np.float32)
S = read_stack(a.stitch + '/avgs.mrcs'); M = np.load(a.stitch + '/avgs_meta.npz'); N = int(M['N']); full, even, odd = S[:N], S[N:2 * N], S[2 * N:3 * N]; pop = M['pop'].astype(float); e1, e2 = M['e1'], M['e2']
n = full.shape[1]; c = a.crop; assert c % 2 == 0 and c <= n
# 2D: soft mask, Fourier crop to c x c; shells on the cropped grid in units of the ORIGINAL box (so lp/hp in A map directly)
x = np.arange(n) - n // 2; rr2 = np.hypot(x[:, None], x[None, :]); R = a.mskdiam / 2 / a.smpd; msk2 = np.clip((R + 4 - rr2) / 8., 0, 1).astype(np.float32)
sl = np.r_[0:c // 2, n - c // 2:n]
def fft2c(img): F = np.fft.fft2(np.fft.ifftshift(img * msk2)); return F[np.ix_(sl, sl)]   # origin-centred transform: a centred object's FFT carries (-1)^(h+k) phases that a tilted slice cannot share
fc = np.fft.fftfreq(c) * c / n; fy2, fx2 = np.meshgrid(fc, fc, indexing='ij'); sh2 = np.round(np.hypot(fx2, fy2) * n).astype(int)   # image arrays are [y, x]; volumes [z, y, x]
klo, khi = int(round(n * a.smpd / a.hp)), int(round(n * a.smpd / a.lp)); band = (sh2 >= klo) & (sh2 <= khi); print('band: shells %d-%d of the %d box (%.0f-%.0f A), %d Fourier samples per slice' % (klo, khi, n, a.hp, a.lp, band.sum()), flush=True)
FE = np.array([fft2c(im) for im in even]); FO = np.array([fft2c(im) for im in odd]); FF = np.array([fft2c(im) for im in full])
# 3D geometry: rotation matrices of the sub-class directions (e3 = 0); four candidate conventions tested on the consensus
def rot(e1d, e2d, e3d=0.):
    a1, a2, a3 = np.radians([e1d, e2d, e3d]); c1, s1, c2, s2, c3, s3 = np.cos(a1), np.sin(a1), np.cos(a2), np.sin(a2), np.cos(a3), np.sin(a3)
    Rz1 = np.array([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]]); Ry = np.array([[c2, 0, s2], [0, 1, 0], [-s2, 0, c2]]); Rz3 = np.array([[c3, -s3, 0], [s3, c3, 0], [0, 0, 1]])
    return Rz1 @ Ry @ Rz3
def slice_index(Rm):   # 3D grid indices of the slice at rotation Rm: image freq (fx, fy, 0) in (x, y, z) -> 3D freq -> array index (z, y, x)
    kk = np.stack([fx2, a.yflip * fy2, np.zeros_like(fx2)], -1).reshape(-1, 3) * n
    p = kk @ Rm.T if a.conv_T else kk @ Rm; q = np.round(p).astype(int)
    ok = (np.abs(q) <= c // 2 - 1).all(1); q = q % c; return q[:, 2] * c * c + q[:, 1] * c + q[:, 0], ok
cons = read_mrc(a.cons).astype(np.float32); x3 = np.arange(n) - n // 2; rr3 = np.sqrt(x3[:, None, None] ** 2 + x3[None, :, None] ** 2 + x3[None, None, :] ** 2); msk3 = np.clip((R + 4 - rr3) / 8., 0, 1)
sl3 = np.r_[0:c // 2, n - c // 2:n]; CF = np.fft.fftn(np.fft.ifftshift(cons * msk3))[np.ix_(sl3, sl3, sl3)].ravel()
def corr_band(A, B, bm=None):
    bm = band if bm is None else bm; a_, b_ = A[bm], B[bm]; return float(np.real(np.vdot(a_, b_)) / (np.linalg.norm(a_) * np.linalg.norm(b_) + 1e-30))
band_lo = (sh2 >= 2) & (sh2 <= khi)   # convention test: from the lowest shells, where a right geometry gives ~0.9 and a wrong one ~0
best = None
for convT in (False, True):
    for sgn in (1., -1.):
        for yf in (1., -1.):
            a.conv_T = convT; a.yflip = yf; cc = []
            for g in range(min(N, 200)):
                idx, ok = slice_index(rot(sgn * e1[g], e2[g])); Fs = np.where(ok, CF[idx], 0).reshape(c, c); cc.append(corr_band(Fs, FF[g], band_lo))
            cc = float(np.mean(cc)); print('  convention transpose=%s e1sign=%+d yflip=%+d: consensus slice vs average corr %.3f' % (convT, int(sgn), int(yf), cc))
            if best is None or cc > best[0]: best = (cc, convT, sgn, yf)
_, a.conv_T, sgn, a.yflip = best; print('chosen: transpose=%s e1sign=%+d yflip=%+d (corr %.3f)' % (a.conv_T, int(sgn), int(a.yflip), best[0]), flush=True)
if best[0] < 0.3: print('GEOMETRY NOT MATCHED (corr < 0.3): the averages are not projections of the consensus at these directions in any tested convention'); sys.exit(2)
if os.environ.get('DIAG_ONLY', '0') == '1': sys.exit(0)
IDX, OK = zip(*[slice_index(rot(sgn * e1[g], e2[g])) for g in range(N)]); IDX, OK = np.array(IDX), np.array(OK)
def insert(F, members, w):   # member-weighted nearest-voxel insertion of slices into a c^3 grid (+ Friedel mates)
    num = np.zeros(c ** 3, np.complex64); den = np.zeros(c ** 3, np.float32)
    for g, wg in zip(members, w):
        v = F[g].ravel() * OK[g]; np.add.at(num, IDX[g], wg * v); np.add.at(den, IDX[g], wg * OK[g])
    V = num / np.maximum(den, 1e-3 * den.max() + 1e-12); Vc = np.conj(V.reshape(c, c, c))[::-1, ::-1, ::-1]; Vc = np.roll(Vc, 1, (0, 1, 2)).ravel()   # Friedel symmetrisation
    return 0.5 * (V + Vc)
def extract(V, g): return np.where(OK[g], V[IDX[g]], 0).reshape(c, c)
# deviations: every average minus the consensus slice at its direction (scale fitted in the band); the shared consensus must not drive the assignment,
# and states are deviation volumes shrunk toward zero where uncovered (consensus + deviation), so coverage cannot be rewarded
CS = np.array([np.where(OK[g], CF[IDX[g]], 0).reshape(c, c) for g in range(N)])
def devs(F):
    out = np.empty_like(F)
    for g in range(N):
        a_, b_ = CS[g][band], F[g][band]; sc = np.real(np.vdot(a_, b_)) / (np.real(np.vdot(a_, a_)) + 1e-30); out[g] = F[g] / max(sc, 1e-12) - CS[g]
    return out
DE, DO = devs(FE), devs(FO); parent = M['parent'].astype(int); bins = np.unique(parent)
if os.environ.get('CENTRE', '1') == '1':   # bin-centred: subtract every bin's mean deviation (CTF, scale, direction error shared by its sub-classes), keep only the within-bin contrasts
    for b in bins:
        mem = np.where(parent == b)[0]; w_ = pop[mem] / pop[mem].sum(); DE[mem] -= (w_[:, None, None] * DE[mem]).sum(0); DO[mem] -= (w_[:, None, None] * DO[mem]).sum(0)
print('deviations from the consensus slices (bin-centred: %s): %d sub-classes in %d bins' % (os.environ.get('CENTRE', '1'), N, len(bins)), flush=True)
from scipy.optimize import linear_sum_assignment
def insert_dev(F, members, w):   # deviation volume: weighted mean of the inserted slices, shrunk to zero where the coverage is thin
    num = np.zeros(c ** 3, np.complex64); den = np.zeros(c ** 3, np.float32)
    for g, wg in zip(members, w):
        v = F[g].ravel() * OK[g]; np.add.at(num, IDX[g], wg * v); np.add.at(den, IDX[g], wg * OK[g])
    eps = 0.05 * np.median(den[den > 0]) if (den > 0).any() else 1.; V = num / (den + eps); Vc = np.conj(V.reshape(c, c, c))[::-1, ::-1, ::-1]; Vc = np.roll(Vc, 1, (0, 1, 2)).ravel()
    return 0.5 * (V + Vc)
K = a.k; cap = int(np.ceil(max(np.bincount(parent).max(), 1) / K)); lab = np.zeros(N, int)
for b in bins: mem = np.where(parent == b)[0]; lab[mem] = rng.permutation(np.arange(len(mem)) % K)   # start: every bin's sub-classes spread over the states
print('E-step: per-bin balanced assignment (at most %d sub-classes of a bin per state), cross-half deviation correlation in the band' % cap, flush=True)
for it in range(a.iters):
    VE = [insert_dev(DE, np.where(lab == k)[0], pop[lab == k]) for k in range(K)]; VO = [insert_dev(DO, np.where(lab == k)[0], pop[lab == k]) for k in range(K)]
    sc = np.zeros((N, K))
    for g in range(N):
        for k in range(K): sc[g, k] = 0.5 * (corr_band(extract(VO[k], g), DE[g]) + corr_band(extract(VE[k], g), DO[g]))   # cross-half: even deviation vs odd volume and vice versa
    new = lab.copy()
    for b in bins:   # balanced assignment inside every bin: rows = its sub-classes, columns = states repeated cap times
        mem = np.where(parent == b)[0]; cost = -np.repeat(sc[mem], cap, axis=1); r_, c_ = linear_sum_assignment(cost); new[mem[r_]] = c_ // cap
    changed = int((new != lab).sum()); lab = new
    print('it %2d: %d sub-classes changed state; populations %s; mean best score %.3f' % (it + 1, changed, ' '.join(str(int(pop[lab == k].sum())) for k in range(K)), sc.max(1).mean()), flush=True)
    if changed == 0 and it > 2: break
# reproducibility of every state (even vs odd volume in the band) and the state maps
f3 = np.fft.fftfreq(c) * c / n; k3 = np.round(np.sqrt(sum(v ** 2 for v in np.meshgrid(f3, f3, f3, indexing='ij'))) * n).astype(int).ravel(); b3 = (k3 >= klo) & (k3 <= khi)
rep = []
for k in range(K):
    e_, o_ = VE[k][b3], VO[k][b3]; rep.append(float(np.real(np.vdot(e_, o_)) / (np.linalg.norm(e_) * np.linalg.norm(o_) + 1e-30)))
    V = (CF + 0.5 * (VE[k] + VO[k])).reshape(c, c, c); write_mrc('%s/state_%02d.mrc' % (a.out, k + 1), np.fft.fftshift(np.real(np.fft.ifftn(V))).astype(np.float32), a.smpd * n / c)
print('states (rep = cross-half correlation of the DEVIATION volume in the band): ' + ' '.join('%d:(%d ptcls, %d sub-classes, rep %.2f)' % (k + 1, int(pop[lab == k].sum()), int((lab == k).sum()), rep[k]) for k in range(K)))
np.savetxt(a.out + '/subclass_state.txt', np.c_[M['gid'], lab + 1], fmt='%d')
L = np.loadtxt(a.clsexp + '/cls_expansion_labels.txt', comments='#').astype(int)   # particle_index parent global_subclass
npart = max(int(L[:, 0].max()), a.npart); pl = np.zeros(npart, int); pl[L[:, 0] - 1] = lab[L[:, 2] - 1] + 1; np.savetxt(a.out + '/particle_states.txt', pl, fmt='%d')
print('particle labels -> %s/particle_states.txt (%d labelled of %d)' % (a.out, int((pl > 0).sum()), npart))

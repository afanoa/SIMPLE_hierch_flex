#!/usr/bin/env python3
"""ManifoldEM-style stitching of direction-bin sub-classes into 3D states by NEIGHBOURING bins only: inside each bin the
sub-class averages are centred (bin mean removed), band-passed and masked; for every pair of neighbouring bins (angular
distance < --maxang) the in-plane rotation between the two bin frames is calibrated on the bin means, then the sub-class
contrasts of the two bins are matched by a cross-half Hungarian assignment (even of A vs odd of B and vice versa).
Labels (K = sub-classes per bin) are propagated from the most populated bin along a maximum-confidence spanning tree.
usage: neighbor_propagate.py STITCH_DIR CLSEXP_DIR BINS_NPZ CAVGS_BINS OUT [--nn 6 --maxang 30 --crop 80 --lp 12 --hp 60]"""
import sys, os, argparse, heapq, numpy as np
from scipy.ndimage import rotate as ndrotate
from scipy.optimize import linear_sum_assignment
ap = argparse.ArgumentParser(); ap.add_argument('stitch'); ap.add_argument('clsexp'); ap.add_argument('bins'); ap.add_argument('cavgs'); ap.add_argument('out')
ap.add_argument('--nn', type=int, default=6); ap.add_argument('--maxang', type=float, default=30.); ap.add_argument('--crop', type=int, default=80); ap.add_argument('--lp', type=float, default=12.); ap.add_argument('--hp', type=float, default=60.)
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--mskdiam', type=float, default=230.); ap.add_argument('--npart', type=int, default=0); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def read_stack(fn):
    h = np.frombuffer(open(fn, 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
    return np.array(np.memmap(fn, np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)), np.float32)
S = read_stack(a.stitch + '/avgs.mrcs'); M = np.load(a.stitch + '/avgs_meta.npz'); N = int(M['N']); full, even, odd = S[:N], S[N:2 * N], S[2 * N:3 * N]; pop = M['pop'].astype(float); parent = M['parent'].astype(int)
B = np.load(a.bins); C = B['centres']; nb = len(C); binavg = read_stack(a.cavgs)
n = full.shape[1]; c = a.crop; sl = np.r_[0:c // 2, n - c // 2:n]
f = np.fft.fftfreq(n); fy, fx = np.meshgrid(f, f, indexing='ij'); k = np.hypot(fx, fy) / a.smpd   # 1/A
bp = (np.exp(-0.5 * (k * a.lp / 1.2) ** 2 * 2.) * (1 - np.exp(-0.5 * (k * a.hp / 1.2) ** 2 * 2.))).astype(np.float32)   # band-pass hp..lp A
def prep(img):   # band-pass, Fourier-crop to c x c, soft mask
    F = np.fft.fft2(np.fft.ifftshift(img)) * bp; F = F[np.ix_(sl, sl)]; r = np.fft.fftshift(np.real(np.fft.ifft2(F))).astype(np.float32); return r
x = np.arange(c) - c // 2; rr = np.hypot(x[:, None], x[None, :]); R = a.mskdiam / 2 / a.smpd * c / n; msk = np.clip((R + 2 - rr) / 4., 0, 1).astype(np.float32)
PF, PE, PO = [np.array([prep(im) for im in X]) for X in (full, even, odd)]; BA = np.array([prep(im) for im in binavg])
members = {b: np.where(parent == b)[0] for b in range(1, nb + 1)}
def par_of(g): return int(parent[g]) - 1
def centred(P):   # remove every bin's pop-weighted mean average
    D = P.copy()
    for b, mem in members.items():
        if len(mem): w = pop[mem] / max(pop[mem].sum(), 1e-12); D[mem] -= (w[:, None, None] * P[mem]).sum(0)
    return D
DF, DE, DO = centred(PF), centred(PE), centred(PO)
def corr(u, v): u = u * msk; v = v * msk; u = u - u.mean(); v = v - v.mean(); return float((u * v).sum() / (np.sqrt((u * u).sum() * (v * v).sum()) + 1e-30))
# neighbour pairs by angular distance between bin centres
ang = np.degrees(np.arccos(np.clip(C @ C.T, -1, 1))); pairs = []
for i in range(nb):
    nnb = [j for j in np.argsort(ang[i])[1:a.nn + 1] if ang[i, j] <= a.maxang]
    pairs += [(i, j) for j in nnb if i < j or i not in [jj for jj in np.argsort(ang[j])[1:a.nn + 1]]]
pairs = sorted(set(tuple(sorted(p)) for p in pairs)); print('%d bins, %d neighbour pairs (nn %d, < %.0f deg), %d sub-classes' % (nb, len(pairs), a.nn, a.maxang, N), flush=True)
angs = np.arange(-180, 180, 4.)
def calib(i, j):   # in-plane rotation taking bin j's frame onto bin i's, from the bin means
    best = (-2, 0.)
    for psi in angs:
        cc = corr(BA[i], ndrotate(BA[j], psi, reshape=False, order=1)); best = max(best, (cc, psi))
    return best[1], best[0]
edges = []
for (i, j) in pairs:
    mi, mj = members[i + 1], members[j + 1]
    if len(mi) < 2 or len(mj) < 2: continue
    psi, cmean = calib(i, j); RJE = np.array([ndrotate(DE[g], psi, reshape=False, order=1) for g in mj]); RJO = np.array([ndrotate(DO[g], psi, reshape=False, order=1) for g in mj])
    Sm = np.zeros((len(mi), len(mj)))
    for u, g in enumerate(mi):
        for v in range(len(mj)): Sm[u, v] = 0.5 * (corr(DE[g], RJO[v]) + corr(DO[g], RJE[v]))
    r_, c_ = linear_sum_assignment(-Sm); conf = float(Sm[r_, c_].mean()); off = float((Sm.sum() - Sm[r_, c_].sum()) / max(Sm.size - len(r_), 1))
    edges.append((conf, i, j, dict(zip(mi[r_], mj[c_])), psi, cmean, off))
conf = np.array([e[0] for e in edges]); print('edges: %d; bin-mean calibration corr median %.2f; matched contrast similarity median %.2f (unmatched %.2f); edges with matched > 0.2: %d' % (len(edges), np.median([e[5] for e in edges]), np.median(conf), np.median([e[6] for e in edges]), int((conf > 0.2).sum())), flush=True)
# propagation along a maximum-confidence spanning tree from the most populated bin
K = max(len(m) for m in members.values()); lab = np.zeros(N, int); seed = max(members, key=lambda b: pop[members[b]].sum()) - 1
for r, g in enumerate(sorted(members[seed + 1], key=lambda g: -pop[g])): lab[g] = r + 1
done = {seed}; adj = {}
for e in edges: adj.setdefault(e[1], []).append(e); adj.setdefault(e[2], []).append(e)
heap = [(-e[0], e[1], e[2], id(e), e) for e in adj.get(seed, [])]; heapq.heapify(heap); nprop = 0
while heap:
    _, i, j, _, e = heapq.heappop(heap); src, dst = (i, j) if i in done and j not in done else ((j, i) if j in done and i not in done else (None, None))
    if src is None: continue
    mp = e[3] if src == i else {v: u for u, v in e[3].items()}   # src sub-class -> dst sub-class
    for gs, gd in mp.items():
        if lab[gs]: lab[gd] = lab[gs]
    for g in members[dst + 1]:   # unmatched sub-classes (bins with fewer members) get a free label
        if lab[g] == 0: lab[g] = 1 + np.argmin(np.bincount(lab[members[dst + 1]], minlength=K + 1)[1:])
    done.add(dst); nprop += 1
    for e2 in adj.get(dst, []): heapq.heappush(heap, (-e2[0], e2[1], e2[2], id(e2), e2))
print('tree propagation: %d of %d bins (K = %d labels); sub-classes per label: %s' % (len(done), nb, K, np.bincount(lab, minlength=K + 1)[1:].tolist()), flush=True)
if os.environ.get('SYNC', '1') == '1':   # permutation synchronisation over ALL edges (spectral, confidence-weighted): robust to wrong edges
    slot = {}
    for b, mem in members.items():
        for r, g in enumerate(sorted(mem, key=lambda g: -pop[g])[:K]): slot[g] = r
    Wm = np.zeros((nb * K, nb * K))
    for e in edges:
        ci, i, j, mp = e[0], e[1], e[2], e[3]; w = max(ci, 0.)
        for gs, gd in mp.items():
            if gs in slot and gd in slot: Wm[i * K + slot[gs], j * K + slot[gd]] += w; Wm[j * K + slot[gd], i * K + slot[gs]] += w
    Wm += np.eye(nb * K); vals, vecs = np.linalg.eigh(Wm); U = vecs[:, -K:]
    Us = U[seed * K:(seed + 1) * K]; lab2 = np.zeros(N, int); nfix = 0
    for b in range(nb):
        Ub = U[b * K:(b + 1) * K]; Cm = Ub @ Us.T; r_, c_ = linear_sum_assignment(-Cm)   # slot r of bin b <-> slot c of the seed
        for g, r in slot.items():
            if par_of(g) == b: lab2[g] = c_[list(r_).index(r)] + 1 if r in r_ else 0
    for g in range(N):
        if lab2[g] == 0: lab2[g] = lab[g]
    print('permutation synchronisation (top-%d eigenvalues %s): sub-classes per label %s' % (K, ' '.join('%.1f' % v for v in vals[-K:]), np.bincount(lab2, minlength=K + 1)[1:].tolist()), flush=True)
    lab = lab2
np.savez(a.out + '/edges.npz', conf=conf, pairs=np.array([(e[1], e[2]) for e in edges]), psi=np.array([e[4] for e in edges]), matches=np.array([(u, v) for e in edges for u, v in e[3].items()]), edge_of_match=np.array([ie for ie, e in enumerate(edges) for _ in e[3]]))
np.savetxt(a.out + '/subclass_state.txt', np.c_[M['gid'], lab], fmt='%d')
L = np.loadtxt(a.clsexp + '/cls_expansion_labels.txt', comments='#').astype(int); npart = max(int(L[:, 0].max()), a.npart); pl = np.zeros(npart, int); pl[L[:, 0] - 1] = lab[np.minimum(L[:, 2] - 1, N - 1)]
np.savetxt(a.out + '/particle_states.txt', pl, fmt='%d'); print('particle labels -> %s/particle_states.txt' % a.out)

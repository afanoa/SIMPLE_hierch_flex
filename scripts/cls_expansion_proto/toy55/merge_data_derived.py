#!/usr/bin/env python3
"""Data-derived merge of leaf states (region-local, spectrally whitened): two states (or pooled groups) are one state unless the DIFFERENCE between them is
itself resolved: the half-map FSC of (dev_even_A - dev_even_B) against (dev_odd_A - dev_odd_B) must reach FLOOR A (the
same floor a leaf needs to count as supported). Greedy: merge the pair whose difference is least resolved, pool their
halves (particle-weighted), repeat until every remaining pair's difference reaches the floor. No similarity cut.
usage: merge_data_derived.py RUNDIR CONSENSUS OUT [--floor 12] ; env LEAFGLOB, FSCMAX, MINREC as merge_states.py"""
import sys, os, glob, re, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc; from merge_states_lib import leaf_ok
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--floor', type=float, default=float(os.environ.get('FLOOR', '12')))
ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--crop', type=int, default=80, help='Fourier box kept (80 = shells to 2.9 A, far beyond any floor)'); a = ap.parse_args()
cons = read_mrc(a.cons).astype(np.float32); n = cons.shape[0]; x = np.arange(n) - n // 2
rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (a.smpd * 10. / 1.2) ** 2 / a.smpd ** 2)
lpc = np.real(np.fft.ifftn(np.fft.fftn(cons) * G)); supp = lpc > 0.2 * lpc.max()
c = a.crop; sl = np.r_[0:c // 2, n - c // 2:n]
Gm = np.exp(-3. * (kx**2 + ky**2 + kz**2) / (a.smpd / 15.) ** 2)   # mask low-pass: 5 % at 15 A (G above is a gentle Gaussian, 60 % at 8 A, useless for masks)
def crop_fft(v):   # masked volume -> Fourier box c^3 (shell indices preserved)
    F = np.fft.fftn(v * msk); return F[np.ix_(sl, sl, sl)]
fc = np.fft.fftfreq(c) * c / n; kx, ky, kz = np.meshgrid(fc, fc, fc, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2).ravel(); nsh = n // 2 + 1
def dev(fn):
    v = read_mrc(fn).astype(np.float32); A = np.stack([lpc[supp], np.ones(supp.sum())], 1); ab = np.linalg.lstsq(A, v[supp], rcond=None)[0]
    return ((v - ab[1]) / max(ab[0], 1e-12) - lpc).astype(np.float32)
names, FE, FO, cnt = [], [], [], []
for d in sorted(glob.glob(a.run + '/fam?/' + os.environ.get('LEAFGLOB', 'region*') + '/level1/maps')):
    m = re.search(r'fam(\d)/(?:zoom)?region([^/]+)/', d); tag = 'f%s_%s' % (m.group(1), m.group(2)); L1 = os.path.dirname(d)   # names as merge_states.py (consolidate_states.py parses them)
    nrec = {}
    if os.path.isfile(L1 + '/clusters.txt'):
        for l in open(L1 + '/clusters.txt'):
            if not l.startswith('#'): nrec[int(l.split()[0])] = int(l.split()[3])
    for fe in sorted(glob.glob(d + '/recvol_state??_even.mrc')):
        s = int(re.findall(r'state(\d+)_even', fe)[0]); fo = fe.replace('_even', '_odd')
        if not os.path.isfile(fo) or not leaf_ok(L1, s): continue
        names.append('%s_s%02d' % (tag, s)); FE.append(crop_fft(dev(fe))); FO.append(crop_fft(dev(fo))); cnt.append(nrec.get(s, 1))
N = len(names); FE, FO, cnt = np.array(FE), np.array(FO), np.array(cnt, float); print('%d leaf states (floor %.0f A)' % (N, a.floor), flush=True)
nsh_cnt = np.bincount(sh, minlength=n // 2 + 1).astype(float)
for i in (range(N) if os.environ.get('WHITEN', '0') == '1' else []):   # optional spectral whitening per state (own mean radial amplitude): a B-factor or scale mismatch between two states is not a difference
    amp = np.sqrt(0.5 * (np.bincount(sh, (np.abs(FE[i]) ** 2).ravel(), n // 2 + 1) + np.bincount(sh, (np.abs(FO[i]) ** 2).ravel(), n // 2 + 1)) / np.maximum(nsh_cnt, 1)); w = (1. / np.maximum(amp, 1e-30))[sh].reshape(FE[i].shape)
    FE[i] *= w; FO[i] *= w
def fsc_of(E, O):
    num = np.bincount(sh, np.real(np.conj(E) * O).ravel(), nsh); d1 = np.bincount(sh, (np.abs(E) ** 2).ravel(), nsh); d2 = np.bincount(sh, (np.abs(O) ** 2).ravel(), nsh)
    return num / np.maximum(np.sqrt(d1 * d2), 1e-30)
def res_of(E, O):   # FSC 0.143 crossing of two half-volumes (Fourier boxes) in A
    fsc = fsc_of(E, O); below = np.where(fsc[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * a.smpd / k
k_hi, k_lo = int(round(n * a.smpd / a.floor)), int(round(n * a.smpd / (1.25 * a.floor)))   # shells of the floor band
def band_fsc(E, O): return float(fsc_of(E, O)[k_lo:k_hi + 1].mean())
def crop_real(v):   # full-box real volume -> real volume on the c^3 grid (Fourier crop)
    F = np.fft.fftn(v)[np.ix_(sl, sl, sl)]; return np.real(np.fft.ifftn(F)).astype(np.float32) * (c / n) ** 3
envc = np.clip(crop_real(np.real(np.fft.ifftn(np.fft.fftn((lpc > 0.1 * lpc.max()).astype(np.float32)) * Gm))), 0, 1); _mc = {}   # soft (10 A) envelope: a sharp mask leaks low-resolution signal into every shell
def region_mask(name):   # the region mask of the board the leaf came from (fam/regions{,_motion}/region_N.mrc), 10 A low-passed, inside the envelope
    m = re.match(r'f(\d)_(motion_)?_?(\d+)', name); fam, mot, reg = m.group(1), m.group(2), m.group(3)
    f = '%s/fam%s/%s/region_%s.mrc' % (a.run, fam, 'regions_motion' if mot else 'regions', reg)
    if f not in _mc:
        if os.path.isfile(f): v = read_mrc(f).astype(np.float32); v = np.real(np.fft.ifftn(np.fft.fftn(v) * Gm)); _mc[f] = np.clip(crop_real(v), 0, 1) * envc
        else: print('  region mask %s MISSING: whole envelope' % f); _mc[f] = envc
    return _mc[f]
RE = np.array([np.real(np.fft.ifftn(F)).astype(np.float32) for F in FE]); RO = np.array([np.real(np.fft.ifftn(F)).astype(np.float32) for F in FO]); del FE, FO
sph = np.clip((0.45 * c - np.sqrt(sum((np.arange(c)[(slice(None),) + (None,) * k + (None,) * (2 - k)] - c // 2) ** 2 for k in range(3))) ) / (0.05 * c), 0, 1).astype(np.float32)
MK = [region_mask(nm) for nm in names] if os.environ.get('MASK', 'region') == 'region' else [sph] * N; print('region masks: %d distinct, mean coverage of the envelope %.2f' % (len(_mc), np.mean([m.sum() / envc.sum() for m in MK])), flush=True)
self_res = np.array([res_of(np.fft.fftn(RE[i] * MK[i]), np.fft.fftn(RO[i] * MK[i])) for i in range(N)]); print('own resolution of the leaves: median %.1f A, worst %.1f A' % (np.median(self_res), self_res.max()))
groups = {i: [i] for i in range(N)}; GE = {i: RE[i].copy() for i in range(N)}; GO = {i: RO[i].copy() for i in range(N)}; GN = {i: cnt[i] for i in range(N)}; GM = {i: MK[i] for i in range(N)}
R = {}
def dres(i, j):   # FSC of the difference between two (pooled) states in the floor band, inside the union of their regions
    m = np.maximum(GM[i], GM[j]); return band_fsc(np.fft.fftn((GE[i] - GE[j]) * m), np.fft.fftn((GO[i] - GO[j]) * m))
def _dres_pair(p): return p, dres(*p)
from multiprocessing import Pool; nproc = int(os.environ.get('NPROC', '16'))
keys = sorted(groups); pairs = [(i, j) for ii, i in enumerate(keys) for j in keys[ii + 1:]]
with Pool(nproc) as pool:
    for p, v in pool.imap_unordered(_dres_pair, pairs, chunksize=64): R[p] = v
for (i, j), v in sorted(R.items(), key=lambda kv: kv[1])[:6]: print('  least-resolved difference at the start: %s vs %s: band FSC %.2f' % (names[i], names[j], v))
curve = lambda E, O: ' '.join('%.2f' % v for v in fsc_of(E, O)[2:k_hi + 6:2])
for pair in (('f1_motion_1_s05_s01', 'f1_motion_1_s05_s03'), ('f1_motion_1_s05', 'f1_motion_1_s01'), ('f1__1_s01', 'f1__1_s02')):
    if pair[0] in names and pair[1] in names:
        i, j = names.index(pair[0]), names.index(pair[1]); m = np.maximum(MK[i], MK[j])
        print('  check %s vs %s: difference band FSC %.2f (floor shells %d-%d); FSC curve of the difference (shells 2..%d step 2): %s' % (pair[0], pair[1], R[tuple(sorted((i, j)))], k_lo, k_hi, k_hi + 5, curve(np.fft.fftn((RE[i] - RE[j]) * m), np.fft.fftn((RO[i] - RO[j]) * m))))
        print('    own FSC curve of %s: %s' % (pair[0], curve(np.fft.fftn(RE[i] * MK[i]), np.fft.fftn(RO[i] * MK[i]))))
print('pairs at the start: difference band FSC median %.2f, %d of %d pairs unresolved (< 0.143) at %.0f A' % (np.median(list(R.values())), sum(v < 0.143 for v in R.values()), len(R), a.floor), flush=True)
if os.environ.get('DIAG_ONLY', '0') == '1': sys.exit(0)
log = []
while True:
    (i, j), worst = min(R.items(), key=lambda kv: kv[1])
    if worst >= 0.143: break
    wi, wj = GN[i], GN[j]; GE[i] = (wi * GE[i] + wj * GE[j]) / (wi + wj); GO[i] = (wi * GO[i] + wj * GO[j]) / (wi + wj); GN[i] = wi + wj; groups[i] += groups.pop(j)
    GM[i] = np.maximum(GM[i], GM[j]); del GE[j], GO[j], GN[j], GM[j]; R = {k: v for k, v in R.items() if j not in k}
    todo = [tuple(sorted((i, k))) for k in groups if k != i]
    with Pool(nproc) as pool:
        for p, v in pool.imap_unordered(_dres_pair, todo, chunksize=8): R[p] = v
    log.append((worst, len(groups)))
lab = np.zeros(N, int)
for g, (k, mem) in enumerate(sorted(groups.items(), key=lambda kv: -GN[kv[0]])):
    for m in mem: lab[m] = g + 1
K = lab.max(); print('%d merges -> %d states; the last merges happened at difference band FSC %s' % (len(log), K, ' '.join('%.2f' % w for w, _ in log[-5:])))
rem = sorted(R.values()); print('remaining pairs: difference band FSC median %.2f, least-resolved pair %.2f (threshold 0.143 at %.0f A)' % (np.median(rem), min(rem), a.floor) if rem else 'one state')
np.savez(a.out + '.npz', names=np.array(names), labels=lab, floor=a.floor, self_res=self_res)
with open(a.out + '.txt', 'w') as o:
    o.write('# data-derived merge: states whose difference half-maps have FSC < 0.143 in the %.0f A band (%s mask, whiten %s) are one state; %d leaves -> %d states\n' % (a.floor, os.environ.get('MASK', 'region'), os.environ.get('WHITEN', '0'), N, K))
    for g in range(1, K + 1):
        mem = np.where(lab == g)[0]; o.write('%d %d %s\n' % (g, int(cnt[mem].sum()), ' '.join(names[m] for m in mem)))
print('-> %s.txt' % a.out)

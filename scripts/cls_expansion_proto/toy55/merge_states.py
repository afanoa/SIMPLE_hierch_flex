#!/usr/bin/env python3
"""Merge tree over all leaf state maps of a run from their even/odd reconstructions: deviation = scale-matched map minus
consensus, low-passed, inside the consensus envelope; similarity(i, j) = mean of corr(dev_even_i, dev_odd_j) and
corr(dev_odd_i, dev_even_j) (cross-half, so noise does not count); average-linkage tree; the state count at each cut and
the merges, written to OUT.txt; the self cross-half reproducibility of every map is in the table.
usage: merge_states.py RUNDIR CONSENSUS OUT [--lp 10]"""
import sys, os, glob, re, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
from scipy.cluster.hierarchy import linkage, fcluster; from scipy.spatial.distance import squareform

def leaf_ok(level1_dir, s, smpd=1.464, fscmax=float(__import__('os').environ.get('FSCMAX', '12')), minrec=int(__import__('os').environ.get('MINREC', '1000'))):
    """leaf support = the state's own half-map FSC (<= fscmax A) and its reconstructed count (>= minrec); junk reproduces across halves but does not reconstruct"""
    import os, numpy as np
    fe, fo = '%s/maps/recvol_state%02d_even.mrc' % (level1_dir, s), '%s/maps/recvol_state%02d_odd.mrc' % (level1_dir, s)
    if not (os.path.isfile(fe) and os.path.isfile(fo)): return True
    nrec = 0
    if os.path.isfile(level1_dir + '/clusters.txt'):
        for l in open(level1_dir + '/clusters.txt'):
            if not l.startswith('#') and int(l.split()[0]) == s: nrec = int(l.split()[3])
    if nrec and nrec < minrec: return False
    e, o = read_mrc(fe).astype(np.float32), read_mrc(fo).astype(np.float32); n = e.shape[0]; x = np.arange(n) - n // 2
    rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
    E, O = np.fft.fftn(e * msk), np.fft.fftn(o * msk); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2)
    num = np.bincount(sh.ravel(), np.real(np.conj(E) * O).ravel()); d1 = np.bincount(sh.ravel(), np.abs(E).ravel() ** 2); d2 = np.bincount(sh.ravel(), np.abs(O).ravel() ** 2)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < 0.143)[0]; sres = (below[0] + 1) if len(below) else n // 2 - 1
    return n * smpd / max(sres, 1) <= fscmax

ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
cons = read_mrc(a.cons).astype(np.float32); n = cons.shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); G = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (a.lp / 1.2)**2)
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G)); cons = lp(cons); supp = cons > 0.2 * cons.max(); env = lp((cons > 0.1 * cons.max()).astype(np.float32)) > 0.3
regmask_on = os.environ.get('REGMASK', '0') == '1'; _mcache = {}
def region_mask(d):   # the focus mask of the region board this leaf (or zoomed leaf) came from, low-passed, inside the envelope
    m = re.search(r'(fam\d)/(?:zoom)?region_?((?:motion_)?\d+)', d); fam, reg = m.group(1), m.group(2)
    f = '%s/%s/%s/region_%s.mrc' % (a.run, fam, 'regions_motion' if reg.startswith('motion_') else 'regions', reg.replace('motion_', ''))
    if f not in _mcache: _mcache[f] = lp(read_mrc(f).astype(np.float32))[env] if os.path.isfile(f) else np.ones(env.sum(), np.float32); print('  region mask %s: %s, mean inside envelope %.3f' % (f, 'found' if os.path.isfile(f) else 'MISSING', float(_mcache[f].mean())))
    return _mcache[f]
def dev(fn, d=None):
    v = lp(read_mrc(fn).astype(np.float32)); A = np.stack([cons[supp], np.ones(supp.sum())], 1); ab = np.linalg.lstsq(A, v[supp], rcond=None)[0]
    x = ((v - ab[1]) / max(ab[0], 1e-12) - cons)[env]; x = x - x.mean()
    if regmask_on and d is not None: x = x * region_mask(d)
    return x
names, E, O = [], [], []
for d in sorted(glob.glob(a.run + '/fam?/' + os.environ.get('LEAFGLOB', 'region*') + '/level1/maps')):
    m = re.search(r'fam(\d)/(?:zoom)?region([^/]+)/', d); tag = 'f%s_%s' % (m.group(1), m.group(2))
    for fe in sorted(glob.glob(d + '/recvol_state??_even.mrc')):
        s = re.findall(r'state(\d+)_even', fe)[0]; fo = fe.replace('_even', '_odd')
        if not os.path.isfile(fo): continue
        if not leaf_ok(os.path.dirname(d), int(s)): print('  leaf %s_s%s dropped (junk: FSC/count)' % (tag, s)); continue
        names.append('%s_s%s' % (tag, s)); E.append(dev(fe, d)); O.append(dev(fo, d))
E, O = np.array(E), np.array(O); N = len(names); print('%d leaf states' % N)
nrm = lambda X: X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-30); En, On = nrm(E), nrm(O)
S = 0.5 * (En @ On.T + On @ En.T); S = 0.5 * (S + S.T); self_rep = np.diag(S).copy()
D = np.clip(1 - S, 0, 2); np.fill_diagonal(D, 0); Z = linkage(squareform(D, checks=False), method='average')
with open(a.out + '.txt', 'w') as o:
    o.write('# cross-half similarity of leaf state deviations (%.0f A, inside the consensus envelope); self reproducibility median %.2f min %.2f\n' % (a.lp, np.median(self_rep), self_rep.min()))
    for t in (0.5, 0.6, 0.7, 0.8, 0.9):
        lab = fcluster(Z, t=1 - t, criterion='distance'); K = lab.max(); o.write('cut %.1f: %d states\n' % (t, K)); print('similarity %.1f: %d states' % (t, K))
    lab = fcluster(Z, t=0.5, criterion='distance'); o.write('# groups at similarity 0.5 (members merge at that cut):\n')
    for c in range(1, lab.max() + 1):
        mem = np.where(lab == c)[0]
        if len(mem) > 1: o.write('  group %d: %s\n' % (c, ' '.join(names[i] for i in mem)))
    o.write('# per state: name self_rep nearest_state nearest_sim\n')
    for i in range(N):
        j = max((k for k in range(N) if k != i), key=lambda k: S[i, k]); o.write('%s %.2f %s %.2f\n' % (names[i], self_rep[i], names[j], S[i, j]))
np.savez(a.out + '.npz', S=S, names=np.array(names), Z=Z)

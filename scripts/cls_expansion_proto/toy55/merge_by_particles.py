#!/usr/bin/env python3
"""Merge tree over the leaf states of a run by PARTICLE overlap instead of map similarity: similarity(i, j) = Jaccard of the
two leaf states' selected particle sets (states_clusters.txt), average linkage; same output format as merge_states.py
(OUT.txt with the counts per cut, OUT.npz with S, names, Z) so consolidate_states.py can cut it. A small cluster (e.g. the
2 % closed population) that shares few particles with a big open state stays its own state even if their maps look alike.
usage: merge_by_particles.py RUNDIR OUT"""
import sys, os, glob, re, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster; from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc

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

run, out = sys.argv[1:3]; names, sets = [], []
for f in sorted(glob.glob(run + '/fam?/' + __import__('os').environ.get('LEAFGLOB', 'region*') + '/level1/states_clusters.txt')):
    m = re.search(r'fam(\d)/(?:zoom)?region([^/]+)/', f); st = np.loadtxt(f, dtype=int)
    for s in range(1, st.max() + 1):
        rows = np.where(st == s)[0]
        if not len(rows): continue
        if not leaf_ok(os.path.dirname(f), s): print('  leaf f%s_%s_s%02d dropped (junk: FSC/count)' % (m.group(1), m.group(2), s)); continue
        names.append('f%s_%s_s%02d' % (m.group(1), m.group(2), s)); sets.append(set(rows.tolist()))
N = len(names); S = np.eye(N)
for i in range(N):
    for j in range(i + 1, N):
        u = len(sets[i] | sets[j]); S[i, j] = S[j, i] = len(sets[i] & sets[j]) / u if u else 0.
D = np.clip(1 - S, 0, 1); np.fill_diagonal(D, 0); Z = linkage(squareform(D, checks=False), method='average')
with open(out + '.txt', 'w') as o:
    o.write('# particle-overlap (Jaccard) similarity of %d leaf states, average linkage\n' % N)
    for t in (0.2, 0.3, 0.5, 0.7):
        K = fcluster(Z, t=1 - t, criterion='distance').max(); o.write('cut %.1f: %d states\n' % (t, K)); print('jaccard %.1f: %d states' % (t, K))
    o.write('# per state: name size nearest_state jaccard\n')
    for i in range(N):
        j = max((k for k in range(N) if k != i), key=lambda k: S[i, k]); o.write('%s %d %s %.2f\n' % (names[i], len(sets[i]), names[j], S[i, j]))
np.savez(out + '.npz', S=S, names=np.array(names), Z=Z); print('%d leaf states ->' % N, out + '.npz')

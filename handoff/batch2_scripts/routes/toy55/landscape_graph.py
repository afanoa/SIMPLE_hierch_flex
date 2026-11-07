#!/usr/bin/env python3
"""Sample a run's state landscape evenly ALONG its main motion, end to end: per family, the states (leaves) are nodes placed by
their maps' deviation from the family parent inside the family's periphery (scale-matched, low-passed, PCA to VAR of the variance);
a diffusion map (Gaussian kernel, alpha = 1 normalisation so the density of leaves does not bend the coordinate) orders them along
the main path; arc length = distance travelled along that order in map space (running mean of the ordered positions); NSEG
segments of equal arc length (ends included) pool their leaves; an interior segment below BUDGET particles joins its neighbour
along the path, the two end segments keep whatever they hold. No truth, no dataset-specific mask.
Writes OUT/states.txt, OUT/groups.txt, OUT/landscape.txt (group family arc_start arc_end arc_mean psi1 psi2 particles sources).
usage: landscape_graph.py CONS PARENT_PREFIX PERI_PREFIX OUT [--nseg 10 --budget 5000 --lp 12 --var 0.9 --smooth 5]"""
import sys, os, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('parent'); ap.add_argument('peri'); ap.add_argument('out')
ap.add_argument('--nseg', type=int, default=10); ap.add_argument('--budget', type=int, default=5000); ap.add_argument('--lp', type=float, default=12.)
ap.add_argument('--var', type=float, default=0.9); ap.add_argument('--smooth', type=int, default=5); ap.add_argument('--smpd', type=float, default=1.464)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
L1 = a.cons + '/fam0/regionall/level1/maps'; st = np.loadtxt(a.cons + '/states.txt', dtype=int); run = os.path.abspath(a.cons)
while run != '/' and not os.path.isdir(run + '/expand'): run = os.path.dirname(run)
famrow = np.zeros(len(st), int)
for F in (1, 2):
    fs = '%s/expand/family%d_subset.txt' % (run, F)
    if os.path.isfile(fs): famrow[np.loadtxt(fs, dtype=int).ravel() - 1] = F
G = {}
for l in open(a.cons + '/groups.txt'):
    if l.startswith('#') or not l.strip(): continue
    w = l.split(); s = int(w[0]); fam = np.bincount(famrow[st == s], minlength=3)[1:]
    if fam.sum() and os.path.isfile('%s/recvol_state%02d.mrc' % (L1, s)): G[s] = (int(w[1]), w[2], int(np.argmax(fam)) + 1)
n = read_mrc(a.parent.replace('{F}', '1') + '.mrc').shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij')
Glp = np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (a.lp / 1.2) ** 2 / a.smpd ** 2)
def lp(v): return np.real(np.fft.ifftn(np.fft.fftn(v) * Glp))
out_lab = np.zeros_like(st); gid = 0; lines = []; land = []
for F in sorted({v[2] for v in G.values()}):
    P = read_mrc(a.parent.replace('{F}', str(F)) + '.mrc').astype(np.float32); M = np.clip(read_mrc(a.peri.replace('{F}', str(F)) + '.mrc'), 0, 1).astype(np.float32)
    Pl = lp(P); supp = Pl > 0.2 * Pl.max(); ss = [s for s in sorted(G) if G[s][2] == F]; X = []
    for s in ss:
        v = lp(read_mrc('%s/recvol_state%02d.mrc' % (L1, s)).astype(np.float32)); A = np.stack([Pl[supp], np.ones(supp.sum())], 1); c = np.linalg.lstsq(A, v[supp], rcond=None)[0]
        X.append((((v - c[1]) / max(c[0], 1e-12) - Pl) * M)[M > 0.05])
    X = np.array(X); w = np.array([G[s][0] for s in ss], float); mu = (w[:, None] * X).sum(0) / w.sum()
    U, S, Vt = np.linalg.svd((X - mu) * np.sqrt(w / w.sum())[:, None], full_matrices=False); ev = S ** 2 / (S ** 2).sum(); k = int(np.searchsorted(np.cumsum(ev), a.var) + 1); Z = (X - mu) @ Vt[:k].T
    D2 = ((Z[:, None, :] - Z[None]) ** 2).sum(-1); nn = np.sort(D2, 1)[:, min(7, len(ss) - 1)]; eps = np.median(nn)
    Kr = np.exp(-D2 / eps); q = Kr.sum(1); K1 = Kr / np.outer(q, q); d = K1.sum(1); Ms = K1 / np.sqrt(np.outer(d, d))   # alpha = 1, symmetric form
    lam, V = np.linalg.eigh(Ms); o = np.argsort(-lam); lam, V = lam[o], V[:, o]; psi = V / V[:, :1]   # right eigenvectors of the Markov matrix
    p1, p2 = psi[:, 1] * lam[1], (psi[:, 2] * lam[2] if len(ss) > 2 else np.zeros(len(ss)))
    order = np.argsort(p1); Zo = Z[order]; h = max(1, a.smooth // 2); Zs = np.array([Zo[max(0, i - h): i + h + 1].mean(0) for i in range(len(ss))])
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Zs, axis=0), axis=1))]; arc_of = np.empty(len(ss)); arc_of[order] = arc; Lt = arc[-1]
    edges = np.linspace(0, Lt, a.nseg + 1); seg = np.clip(np.searchsorted(edges, arc_of, side='right') - 1, 0, a.nseg - 1)
    segs = [list(np.where(seg == j)[0]) for j in range(a.nseg)]; segs = [g for g in segs if g]
    i = 1   # interior segments below the budget join the neighbour along the path with fewer particles; the two ends stay
    while i < len(segs) - 1:
        if w[segs[i]].sum() < a.budget:
            jn = i - 1 if (i - 1 > 0 and (i + 1 == len(segs) - 1 or w[segs[i - 1]].sum() <= w[segs[i + 1]].sum())) or i + 1 == len(segs) - 1 else i + 1
            if jn in (0, len(segs) - 1) and len(segs) > 3: jn = i + 1 if jn == i - 1 else i - 1
            segs[jn] = segs[jn] + segs[i]; segs.pop(i); i = max(1, i - 1)
        else: i += 1
    print('family %d: %d states, %d components for %.0f %% of the variance; diffusion eigenvalues %s; path length %.1f (map units), %d segments: %s' % (
        F, len(ss), k, 100 * a.var, ' '.join('%.3f' % v for v in lam[1:4]), Lt, len(segs), ' '.join('%d' % int(w[g].sum()) for g in segs)))
    for g in segs:
        gid += 1; mem = [ss[i] for i in g]
        for s in mem: out_lab[st == s] = gid
        lines.append('%d %d %s' % (gid, int(w[g].sum()), '+'.join(G[s][1] for s in mem)))
        land.append('%d %d %.2f %.2f %.2f %.4f %.4f %d %s' % (gid, F, arc_of[g].min(), arc_of[g].max(), (w[g] * arc_of[g]).sum() / w[g].sum(), (w[g] * p1[g]).sum() / w[g].sum(), (w[g] * p2[g]).sum() / w[g].sum(), int(w[g].sum()), ','.join(str(s) for s in mem)))
np.savetxt(a.out + '/states.txt', out_lab, fmt='%d')
with open(a.out + '/groups.txt', 'w') as fo: fo.write('# path state: particles, members (source states); even arc length, ends kept\n'); fo.write('\n'.join(lines) + '\n')
with open(a.out + '/landscape.txt', 'w') as fo: fo.write('# group family arc_start arc_end arc_mean psi1 psi2 particles source_states\n'); fo.write('\n'.join(land) + '\n')
print('%d path states from %d source states' % (gid, len(G)))

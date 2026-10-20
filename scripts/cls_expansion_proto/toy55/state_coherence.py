#!/usr/bin/env python3
"""Are the merged states contiguous regions of the conformer landscape? For a dense landscape (55 conformers, few
particles each) no state can be pure, so purity is the wrong readout. Per merged state: the conformers it holds (share
>= --minshare), the mean pairwise similarity of their docked maps' deviations from the consensus (10 A, inside the
envelope), and the same for random conformer sets of the same size (z-score). Coherent states hold neighbouring
conformers. usage: state_coherence.py CONS_DIR GT_PARTICLES DOCKED_DIR CONSENSUS [--minshare 0.05]"""
import sys, os, glob, re, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
ap = argparse.ArgumentParser(); ap.add_argument('cons'); ap.add_argument('gt'); ap.add_argument('docked'); ap.add_argument('consensus'); ap.add_argument('--minshare', type=float, default=0.05); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args()
C = read_mrc(a.consensus).astype(np.float32); n = C.shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (a.smpd / 10.) ** 2))).astype(np.float32); Cl = lp(C); env = lp((Cl > 0.1 * Cl.max()).astype(np.float32)) > 0.3; supp = Cl > 0.2 * Cl.max()
dev = {}
for fn in sorted(glob.glob(a.docked + '/conf??_*_m160.mrc')):
    c = int(re.search(r'conf(\d+)_', fn).group(1)); v = lp(read_mrc(fn).astype(np.float32)); A = np.stack([Cl[supp], np.ones(supp.sum())], 1); ab = np.linalg.lstsq(A, v[supp], rcond=None)[0]
    d = ((v - ab[1]) / max(ab[0], 1e-12) - Cl)[env]; dev[c] = (d - d.mean()) / (np.linalg.norm(d - d.mean()) + 1e-30)
ids = sorted(dev); D = np.array([dev[c] for c in ids]); S = D @ D.T; idx = {c: i for i, c in enumerate(ids)}
st = np.loadtxt(a.cons + '/states.txt', dtype=int); conf = {int(w[0]): int(w[1]) for w in (l.split() for l in open(a.gt)) if w and not w[0].startswith('#')}; g = np.array([conf.get(i + 1, 0) for i in range(len(st))])
rng = np.random.default_rng(0); K = st.max(); out = []
def mean_pair(sel): return float(np.mean([S[i, j] for ii, i in enumerate(sel) for j in sel[ii + 1:]])) if len(sel) > 1 else float('nan')
for s in range(1, K + 1):
    m = (st == s) & (g > 0); nn = int(m.sum())
    if nn == 0: continue
    cnt = np.bincount(g[m], minlength=max(ids) + 1); share = cnt / nn; members = [c for c in ids if share[c] >= a.minshare]
    if len(members) < 2: out.append((s, nn, len(members), float('nan'), float('nan'))); continue
    coh = mean_pair([idx[c] for c in members]); null = np.array([mean_pair(list(rng.choice(len(ids), len(members), replace=False))) for _ in range(200)]); z = (coh - null.mean()) / (null.std() + 1e-12)
    out.append((s, nn, len(members), coh, z))
print('state  particles  conformers(>=%.0f%%)  mean pairwise conformer similarity  z vs random sets' % (100 * a.minshare))
for s, nn, k, coh, z in sorted(out, key=lambda r: -r[1]): print('%3d  %7d  %3d  %6.3f  %6.2f' % (s, nn, k, coh, z))
zs = np.array([r[4] for r in out if np.isfinite(r[4])]); ks = np.array([r[2] for r in out])
print('summary: %d states; conformers per state median %d; coherence z median %.2f (fraction of states with z > 2: %.2f); random-set baseline similarity %.3f' % (len(out), int(np.median(ks)), np.median(zs) if len(zs) else float('nan'), float(np.mean(zs > 2)) if len(zs) else float('nan'), float(S[np.triu_indices(len(ids), 1)].mean())))

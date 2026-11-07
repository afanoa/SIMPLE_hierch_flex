#!/usr/bin/env python3
"""What is the difference between two merged groups of a split test? D = (small group map) - (large group map) in each
tree (count-weighted candidate maps, each scaled to the consensus), then per shell: FSC(D_A, D_B), |D|^2 / |consensus|^2,
FSC of the even/odd halves of D within a tree; and a PNG of the central slices of D_A, D_B and the consensus (15 A).
usage: diff_diag.py RUN_A RUN_B TAG SMALL_A SMALL_B CONSENSUS OUT_PNG   (SMALL_*: comma-separated candidate state ids)"""
import sys, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc
ra, rb, tag, sa, sb, consf, png = sys.argv[1:8]; smpd = 1.464
cons = read_mrc(consf).astype(np.float64); n = cons.shape[0]; x = np.arange(n) - n // 2
rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2); nsh = n // 2 + 1
G15 = np.exp(-0.5 * (kx**2 + ky**2 + kz**2) * (15. / smpd) ** 2 / 4.); supp = np.real(np.fft.ifftn(np.fft.fftn(cons) * G15)) > 0.2 * cons.max()
def scaled(v): A = np.stack([cons[supp], np.ones(supp.sum())], 1); ab = np.linalg.lstsq(A, v[supp], rcond=None)[0]; return (v - ab[1]) / ab[0]
def group_maps(run, small):
    D = run + '/fam0/region' + tag + '/level1'; st = np.loadtxt(D + '/states.txt').astype(int); K = st.max(); cnt = np.bincount(st, minlength=K + 1)
    small = [int(s) for s in small.split(',')]; out = {}
    for h in ('', '_even', '_odd'):
        M = {s: scaled(read_mrc('%s/maps/recvol_state%02d%s.mrc' % (D, s, h)).astype(np.float64)) for s in range(1, K + 1)}
        S = sum(cnt[s] * M[s] for s in small) / sum(cnt[s] for s in small); L = sum(cnt[s] * M[s] for s in range(1, K + 1) if s not in small) / sum(cnt[s] for s in range(1, K + 1) if s not in small)
        out[h] = (S - L) * msk
    return out, sum(cnt[s] for s in small), cnt[1:].sum() - sum(cnt[s] for s in small)
def shells(a, b=None):
    A = np.fft.fftn(a); B = A if b is None else np.fft.fftn(b)
    num = np.bincount(sh.ravel(), np.real(A * np.conj(B)).ravel(), nsh); da = np.bincount(sh.ravel(), (np.abs(A) ** 2).ravel(), nsh); db = np.bincount(sh.ravel(), (np.abs(B) ** 2).ravel(), nsh)
    return num / np.maximum(np.sqrt(da * db), 1e-30), da
DA, nsA, nlA = group_maps(ra, sa); DB, nsB, nlB = group_maps(rb, sb)
fx, pA = shells(DA['']); _, pB = shells(DB['']); fsc_ab, _ = shells(DA[''], DB['']); fe_a, _ = shells(DA['_even'], DA['_odd']); fe_b, _ = shells(DB['_even'], DB['_odd']); _, pc = shells(cons * msk)
print('groups: A small %d / large %d, B small %d / large %d' % (nsA, nlA, nsB, nlB))
print('%5s %7s %9s %9s %9s %11s %11s' % ('shell', 'res A', 'FSC(DA,DB)', 'eo FSC A', 'eo FSC B', '|DA|2/|C|2', '|DB|2/|C|2'))
for k in list(range(1, 12)) + list(range(12, 45, 4)):
    print('%5d %7.1f %9.3f %9.3f %9.3f %11.2e %11.2e' % (k, n * smpd / k, fsc_ab[k], fe_a[k], fe_b[k], pA[k] / max(pc[k], 1e-30), pB[k] / max(pc[k], 1e-30)))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * G15)); vs = [lp(DA['']), lp(DB['']), lp(cons * msk)]
fig, ax = plt.subplots(3, 3, figsize=(9, 9))
for r, (v, nm) in enumerate(zip(vs, ('D tree A (small - large)', 'D tree B', 'consensus'))):
    c = n // 2; sl = [v[c], v[:, c], v[:, :, c]]; m = np.abs(np.r_[sl[0].ravel(), sl[1].ravel(), sl[2].ravel()]).max()
    for j in range(3): ax[r, j].imshow(sl[j], cmap='RdBu_r' if r < 2 else 'gray', vmin=-m if r < 2 else None, vmax=m if r < 2 else None); ax[r, j].set_xticks([]); ax[r, j].set_yticks([])
    ax[r, 0].set_ylabel(nm)
plt.tight_layout(); plt.savefig(png, dpi=80); print('png', png)

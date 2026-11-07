#!/usr/bin/env python3
"""What separates the states of a labelling besides structure: per state (and per conformer inside it) the refined pose error against
the true poses (projection-direction angle; share > 20 deg), the alignment score (corr), defocus, shift magnitude, and the image
contrast (std inside the mask disc of every image). Rows: project rows = stack rows.
usage: state_nuisance.py LABELS REFINED_FIELDS TRUE_ORITAB STACK GT [--mskdiam 220 --smpd 1.464 --minn 200]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('labels'); ap.add_argument('fields'); ap.add_argument('truth'); ap.add_argument('stack'); ap.add_argument('gt')
ap.add_argument('--mskdiam', type=float, default=220.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--minn', type=int, default=200); a = ap.parse_args()
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
def tab(f, keys): return np.array([[float(kv(l).get(k, 'nan')) for k in keys] for l in open(f) if l.strip()])
K = ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'corr'); R = tab(a.fields, K); T = tab(a.truth, ('e1', 'e2', 'e3', 'x', 'y'))
lab = np.loadtxt(a.labels, dtype=int); G = np.loadtxt(a.gt, comments='#'); conf = np.full(len(lab), -1); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
def dirv(e):
    ph, th = np.radians(e[:, 0]), np.radians(e[:, 1]); return np.c_[np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]
err = np.degrees(np.arccos(np.clip((dirv(R[:, :2]) * dirv(T[:, :2])).sum(1), -1, 1)))
h = np.frombuffer(open(a.stack, 'rb').read(1024), np.int32); nx, nz, nsym = int(h[0]), int(h[2]), int(h[23])
I = np.memmap(a.stack, np.float32, 'r', offset=1024 + nsym, shape=(nz, nx, nx)); x = np.arange(nx) - nx // 2; disc = np.hypot(x[:, None], x[None, :]) < a.mskdiam / 2 / a.smpd
sd = np.zeros(nz, np.float32)
for b in range(0, nz, 5000): blk = np.asarray(I[b:b + 5000], np.float32); sd[b:b + len(blk)] = blk[:, disc].std(1)
i_ = np.arange(40) + 0.5; z_ = 1 - 2 * i_ / 40; r_ = np.sqrt(np.maximum(1 - z_ * z_, 0)); g_ = np.pi * (3 - np.sqrt(5)) * i_
vb = np.argmax(dirv(R[:, :2]) @ np.c_[r_ * np.cos(g_), r_ * np.sin(g_), z_].T, 1); df_ = 0.5 * (R[:, 5] + R[:, 6]); dbin = np.clip(((df_ - np.nanmin(df_)) / (np.nanmax(df_) - np.nanmin(df_) + 1e-9) * 10).astype(int), 0, 9)
allm = lab > 0; hv = np.bincount(vb[allm], minlength=40) / allm.sum(); hd = np.bincount(dbin[allm], minlength=10) / allm.sum()
print('state   n   conformers (top)          pose err med / >20deg   corr     defocus  |shift|  image sd   view TVD  defocus TVD')
for s in np.unique(lab[lab > 0]):
    m = lab == s
    if m.sum() < a.minn: continue
    bc = np.bincount(conf[m][conf[m] >= 0]); top = ' '.join('%d:%.0f%%' % (k, 100. * bc[k] / m.sum()) for k in np.argsort(-bc)[:2] if bc[k] > 0)
    tv = 0.5 * np.abs(np.bincount(vb[m], minlength=40) / m.sum() - hv).sum(); td = 0.5 * np.abs(np.bincount(dbin[m], minlength=10) / m.sum() - hd).sum()
    print('%5d %6d  %-24s  %6.1f / %4.0f%%        %.3f    %.3f    %.2f    %.5f    %.3f     %.3f' % (s, m.sum(), top, np.median(err[m]), 100 * np.mean(err[m] > 20), np.nanmean(R[m, 7]), np.nanmean(0.5 * (R[m, 5] + R[m, 6])), np.mean(np.hypot(R[m, 3], R[m, 4])), np.mean(sd[m]), tv, td))

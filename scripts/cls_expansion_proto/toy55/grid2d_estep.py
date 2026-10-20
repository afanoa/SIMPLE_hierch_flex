#!/usr/bin/env python3
"""2-D deformation latent per particle: references = canonical window map warped along flow modes 1 AND 2 on an
n1 x n2 grid (from deform_proto's flows_dec2.npy / flow_scores.npy), cross-half E-step of a particle subset against
the grid, readouts vs truth: Hungarian accuracy over the grid cells, purity, Spearman of the mode-1 posterior mean
with the truth axis, and the same for the 1-D grid along mode 1 alone (same subset, same band).
usage: grid2d_estep.py DEFORMDIR MAPDIR OUT --nw 8 --info windows.txt --stk --table --gt --states --axis [--n1 7 --n2 5 --lp 12 --nsub 30000]"""
import sys, os, argparse, types, time
import numpy as np
os.environ.setdefault('TOY_SMPD', '1.464'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model_class import read_mrc, write_mrc, warp
from latent_windows import cmd_estep, parse_axis, read_gt
from scipy.ndimage import zoom
from scipy.stats import spearmanr
SMPD = float(os.environ['TOY_SMPD'])

ap = argparse.ArgumentParser(); ap.add_argument('deformdir'); ap.add_argument('mapdir'); ap.add_argument('out'); ap.add_argument('--nw', type=int, default=8); ap.add_argument('--info', required=True)
ap.add_argument('--stk', required=True); ap.add_argument('--table', required=True); ap.add_argument('--gt', required=True); ap.add_argument('--states', required=True); ap.add_argument('--axis', required=True)
ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=5); ap.add_argument('--lp', type=float, default=12.); ap.add_argument('--nsub', type=int, default=30000); ap.add_argument('--mskdiam', type=float, default=200.); ap.add_argument('--tag', default='grid2d')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time(); T = a.tag
Ud = np.load(a.deformdir + '/flows_dec2.npy').astype(float); scores = np.load(a.deformdir + '/flow_scores.npy'); nw = Ud.shape[0]; B2 = Ud.shape[2]; B = 2 * B2
info = np.loadtxt(a.info, comments='#', ndmin=2); n = info[:, 1]; can = a.nw // 2 + 1
X = Ud.reshape(nw, -1); Xc = X - X.mean(0); ubar = Ud.mean(0)
modes = []
for q in range(2):
    sq = scores[:, q]; m = (Xc.T @ sq) / max(np.sum(sq ** 2), 1e-30); modes.append(np.array([zoom(c, 2, order=1) for c in m.reshape(3, B2, B2, B2)]) * 2)
ubar_full = np.array([zoom(c, 2, order=1) for c in ubar]) * 2
def grid_of(sc, k, pad=0.15):
    lo, hi = sc.min(), sc.max(); return np.linspace(lo - pad * (hi - lo), hi + pad * (hi - lo), k)
g1, g2 = grid_of(scores[:, 0], a.n1), grid_of(scores[:, 1], a.n2)
print('%s: %d windows, canonical %d, mode-1 grid %s, mode-2 grid %s (window scores: mode 1 %s, mode 2 %s)' % (T, nw, can, np.array2string(g1, precision=2), np.array2string(g2, precision=2), np.array2string(scores[:, 0], precision=2), np.array2string(scores[:, 1], precision=2)), flush=True)
fn = lambda w, suf='': '%s/recvol_state%02d%s.mrc' % (a.mapdir, w, suf)
Mc = {suf: read_mrc(fn(can, suf)) for suf in ('_even', '_odd')}; Mc[''] = 0.5 * (Mc['_even'] + Mc['_odd'])
# 2-D grid references (cell index = i1 * n2 + i2 + 1) and the 1-D grid along mode 1 (n1 refs, mode 2 at its mean)
for name, cells in (('g2d', [(z1, z2) for z1 in g1 for z2 in g2]), ('g1d', [(z1, scores[:, 1].mean()) for z1 in g1])):
    d = '%s/%s' % (a.out, name); os.makedirs(d, exist_ok=True)
    for k, (z1, z2) in enumerate(cells):
        u = ubar_full + z1 * modes[0] + z2 * modes[1]
        for suf in ('', '_even', '_odd'): write_mrc('%s/recvol_state%02d%s.mrc' % (d, k + 1, suf), warp(Mc[suf], u), SMPD)
    print('%s %s: %d references written (%.0f s)' % (T, name, len(cells), time.time() - t0), flush=True)
base = dict(stk=a.stk, table=a.table, gt=a.gt, refs='', lp=a.lp, states=a.states, axis=a.axis, prev='', prior_sigma=1.5, nsub=a.nsub, mskdiam=a.mskdiam)
cmd_estep(types.SimpleNamespace(mapdir=a.out + '/g1d', nw=a.n1, out=a.out + '/coord_g1d.txt', tag=T + ' 1-D grid', **base))
cmd_estep(types.SimpleNamespace(mapdir=a.out + '/g2d', nw=a.n1 * a.n2, out=a.out + '/coord_g2d.txt', tag=T + ' 2-D grid', **base))
# 2-D readout: posterior over cells -> posterior-mean mode-1 and mode-2 coordinates vs the truth axis / conformer
P = np.load(a.out + '/coord_g2d.txt.post.npz'); post, pinds = P['post'], P['pinds']; gt = read_gt(a.gt); axis = parse_axis(a)
conf = np.array([gt.get(int(p), (0, 0))[0] for p in pinds]); good = conf > 0; ax = np.array([axis.get(c, 0.) for c in conf])
z1 = post @ np.repeat(g1, a.n2); z2 = post @ np.tile(g2, a.n1)
print('%s 2-D posterior means: Spearman(z1, truth axis) %+.3f, Spearman(z2, truth axis) %+.3f; per-conformer mean (z1, z2): %s' %
      (T, spearmanr(z1[good], ax[good])[0], spearmanr(z2[good], ax[good])[0], ' '.join('%d:(%+.2f,%+.2f)' % (c, z1[good][conf[good] == c].mean(), z2[good][conf[good] == c].mean()) for c in np.unique(conf[good]))))
# conformer separability of the 2-D posterior-mean coordinates: pairwise d' in (z1, z2)
Z = np.c_[z1, z2][good]; cg = conf[good]; cls = np.unique(cg); mu = np.array([Z[cg == c].mean(0) for c in cls]); W = sum(np.cov(Z[cg == c].T) * np.sum(cg == c) for c in cls) / len(cg); Wi = np.linalg.inv(W + 1e-9 * np.eye(2))
dp = [np.sqrt((mu[i] - mu[j]) @ Wi @ (mu[i] - mu[j])) for i in range(len(cls)) for j in range(i + 1, len(cls))]
print("%s 2-D coordinate: pairwise conformer d' median %.2f min %.2f max %.2f (1-D along z1 only: median %.2f)" % (T, np.median(dp), np.min(dp), np.max(dp),
      np.median([abs(mu[i, 0] - mu[j, 0]) / np.sqrt(W[0, 0]) for i in range(len(cls)) for j in range(i + 1, len(cls))])))
print('%s DONE (%.0f s)' % (T, time.time() - t0))

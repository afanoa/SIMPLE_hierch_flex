#!/usr/bin/env python3
"""One expansion sub-class against our 3D maps: the sub-class average next to pose-mean renders (projections at the
members' poses, class frame: shift sign -1, scale 1, conv 0, as closed_subclasses.py) of a list of maps, at two bands;
normalised cross-correlation of each render with the average over the whole image and inside the projected closed-arm
region; member statistics against the per-particle focal scores. usage: render_subclass.py SUB META CAVGS OUT [--maps name=path ...]"""
import sys, glob, argparse, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool')
from scipy.ndimage import map_coordinates, gaussian_filter
from view_tangent_toy import euler2m
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('sub', type=int); ap.add_argument('meta'); ap.add_argument('cavgs'); ap.add_argument('out')
ap.add_argument('--maps', nargs='+', default=[]); ap.add_argument('--poses', default='', help='npz with refined e1,e2,e3,x,y per row (used where available)'); ap.add_argument('--scores', default='', help='per-particle focal score file (row score)')
ap.add_argument('--top', default='', help='states file(s) comma-separated: fraction of members with state 1'); ap.add_argument('--nmem', type=int, default=400); a = ap.parse_args()
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; box = 160; rng = np.random.default_rng(0)
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
fftc = lambda v: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v)))
M = np.load(a.meta); flags = M['flags'].astype(int); ref = np.where(flags == 1)[0]; lab = M['lab692'].astype(int); N = len(flags)
E = {k: M[k].astype(float) for k in ('e1', 'e2', 'e3', 'x', 'y')}
if a.poses:
    P = np.load(a.poses)
    if all(k in P for k in ('e1', 'e2', 'e3', 'x', 'y')) and len(P['e1']) == N: E = {k: P[k].astype(float) for k in E}; print('poses: refined table %s' % a.poses)
rows = ref[np.where(lab == a.sub)[0]]; print('sub-class %d: %d members (rows of the 932k project)' % (a.sub, len(rows)))
if a.scores:
    S = np.loadtxt(a.scores, comments='#'); sc = np.full(N + 1, np.nan); sc[S[:, 0].astype(int)] = S[:, 1]; v = sc[rows + 1]; ok = ~np.isnan(v)
    print('  on the family-1 board: %d; focal score median %+.2f, quartiles %+.2f / %+.2f; above +0.84 (the docked top-30k cut): %.0f %%' % (ok.sum(), np.median(v[ok]), np.percentile(v[ok], 25), np.percentile(v[ok], 75), 100 * np.mean(v[ok] > 0.84)))
for sf in [s for s in a.top.split(',') if s]:
    st = np.loadtxt(sf, dtype=int); print('  members in state 1 of %s: %d of %d (%.0f %%)' % (sf.split('/')[-3] + '/' + sf.split('/')[-1], int((st[rows] == 1).sum()), len(rows), 100. * np.mean(st[rows] == 1)))
mem = rows if len(rows) <= a.nmem else rng.choice(rows, a.nmem, replace=False)
Rs = np.array([euler2m([-E['e1'][i], E['e2'][i], -E['e3'][i]]) for i in mem]); xs, ys = E['x'][mem], E['y'][mem]
zz, yy, xx = np.indices((box, box, box)); r3 = np.hypot(np.hypot(xx - box / 2, yy - box / 2), zz - box / 2); msk = 1. / (1. + np.exp((r3 - 116 / SMPD) / 2.))
ky2, kx2 = np.indices((box, box)); kx2 = kx2 - box // 2; ky2 = ky2 - box // 2; r2 = np.hypot(kx2, ky2)
def render(VF, kmax):
    sel = (r2 <= kmax) & (r2 >= 1); kx, ky = kx2[sel].astype(float), ky2[sel].astype(float); c = box // 2; hk = np.stack([kx, ky, np.zeros_like(kx)], 1)
    acc = np.zeros(len(kx), complex)
    for m in range(len(Rs)):
        loc = hk @ Rs[m]; coords = [loc[:, 2] + c, loc[:, 1] + c, loc[:, 0] + c]
        sl = map_coordinates(VF.real, coords, order=1, mode='nearest') + 1j * map_coordinates(VF.imag, coords, order=1, mode='nearest')
        acc += sl * np.exp(-2j * np.pi * (kx * (-xs[m]) + ky * (-ys[m])) / box)
    F2 = np.zeros((box, box), complex); F2[sel] = acc / len(Rs); return np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F2))))
def bandlimit(img, kmax):
    F2 = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img))); F2[(r2 > kmax) | (r2 < 1)] = 0; return np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F2))))
h = np.fromfile(a.cavgs, np.int32, 256); nimg = int(h[2]); stk = np.memmap(a.cavgs, np.float32, 'r', offset=1024 + int(h[23]), shape=(nimg, box, box)); avg = np.array(stk[a.sub - 1], float)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2 * L18.std(); diff = L18 - Lc
arm = ((diff > diff.mean() + 3 * diff.std()) & dens18).astype(float)
maps = [('average', None)] + [(s.split('=')[0], s.split('=', 1)[1]) for s in a.maps]
VFs = {nm: fftc(read_mrc(p) * msk) for nm, p in maps if p}; VFarm = fftc(gaussian_filter(arm, 1.5))
fig, ax = plt.subplots(2, len(maps), figsize=(2.6 * len(maps), 5.6)); out = {}
for row, res in enumerate((8., 12.)):
    kmax = int(np.floor(box * SMPD / res)); A = bandlimit(avg, kmax); armp = render(VFarm, kmax); reg = armp > 0.3 * armp.max()
    An = (A - A.mean()) / A.std()
    for col, (nm, p) in enumerate(maps):
        img = An if p is None else render(VFs[nm], kmax); imn = (img - img.mean()) / max(img.std(), 1e-30)
        cc = float(np.mean(An * imn)); ccr = float(np.corrcoef(An[reg], imn[reg])[0, 1]) if reg.sum() > 10 else float('nan')
        axx = ax[row, col]; axx.imshow(imn, cmap='gray', vmin=-2.5, vmax=4); axx.contour(reg, levels=[0.5], colors='r', linewidths=0.5); axx.axis('off')
        axx.set_title('%s @%g A' % (nm, res) if p is None else '%s\nncc %.2f, arm region %.2f' % (nm, cc, ccr), fontsize=8)
        if p is not None: out.setdefault(nm, []).append((res, cc, ccr))
plt.suptitle('sub-class %d: average vs pose-mean renders at its %d member poses (red: projected closed-arm region)' % (a.sub, len(mem))); plt.tight_layout(); plt.savefig(a.out, dpi=120)
for nm, v in out.items(): print('  %-28s ' % nm + ' | '.join('@%g A: ncc %.2f, arm region %.2f' % t for t in v))
print('figure -> ' + a.out)

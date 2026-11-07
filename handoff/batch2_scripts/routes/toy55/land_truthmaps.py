#!/usr/bin/env python3
"""land_truthmaps.py: ORACLE references for the landscape tests. Truth maps soft-masked the way simulate_particles masks a volume
(mask3D_soft: cosine edge from r0 to box/2 with box = min(nint(2 (msk + 6)), n), msk = mskdiam / 2 / smpd px), plus the lsamp_score2.py
inputs OUT/leaves.txt (name map map count) and OUT/labels.txt (master row -> 1..K = true conformer index, 0 = not scored).
usage: land_truthmaps.py MAPDIR GT OUT [--ids 200:230 --mskdiam 200 --smpd 1.464]"""
import sys, os, glob, argparse, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, write_mrc
ap = argparse.ArgumentParser(); ap.add_argument('mapdir'); ap.add_argument('gt'); ap.add_argument('out'); ap.add_argument('--ids', default='200:230')
ap.add_argument('--mskdiam', type=float, default=200.); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
i0, i1 = [int(x) for x in a.ids.split(':')]; ids = list(range(i0, i1 + 1))
conf = np.loadtxt(a.gt, comments='#')[:, 1].astype(int); lab = np.zeros(len(conf), int); m = (conf >= i0) & (conf <= i1); lab[m] = conf[m] - i0 + 1
n = read_mrc(glob.glob('%s/conf%d_*.mrc' % (a.mapdir, ids[0]))[0]).shape[0]; msk = a.mskdiam / 2. / a.smpd
minlen = min(int(np.rint(2. * (msk + 6.))), n); maxrad = float(minlen // 2); width = 2. * (maxrad - msk); r0 = maxrad - width
x = np.arange(n) - n / 2.; r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
w = np.where(r <= r0, 1., np.where(r >= maxrad, 0., 0.5 * (np.cos((r - r0) / width * np.pi) + 1.)))
L = []
for k in ids:
    fo = '%s/conf%d_m%d.mrc' % (a.out, k, int(a.mskdiam)); write_mrc(fo, read_mrc(glob.glob('%s/conf%d_*.mrc' % (a.mapdir, k))[0]) * w, a.smpd)
    L.append('t%02d %s %s %d' % (k - i0, fo, fo, int((conf == k).sum())))
open(a.out + '/leaves.txt', 'w').write('\n'.join(L) + '\n'); np.savetxt(a.out + '/labels.txt', lab, fmt='%d')
print('%d truth maps masked at mskdiam %.0f (cosine ramp %.1f .. %.1f px); %d of %d rows scored' % (len(ids), a.mskdiam, r0, maxrad, m.sum(), len(conf)))

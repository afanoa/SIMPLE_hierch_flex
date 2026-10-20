#!/usr/bin/env python3
"""Concatenate the per-class simulations in a random order: all.mrcs, all.txt (truth: poses + CTF, state=conformer id,
0 for junk), imp.txt (CTF only), gt_particles.txt (particle conformer closed e1 e2). usage: merge55.py SIMDIR CSV CLOSED_IDS
Streams the images in output chunks (memory ~ CHUNK images), so sets of millions of particles need no RAM for the stack."""
import sys, csv, os, numpy as np
T, csvf, closed = sys.argv[1], sys.argv[2], set(int(v) for v in sys.argv[3].split(','))
rows = list(csv.DictReader(open(csvf))); names = ['conf%02d' % int(r['conformer']) for r in rows] + ['blob', 'broken']; labs = [int(r['conformer']) for r in rows] + [0, 0]
def kv(line): return dict(t.split('=', 1) for t in line.split() if '=' in t)
mm, oris, lab, src = [], [], [], []
for k, (name, l) in enumerate(zip(names, labs)):
    fn = '%s/stk_%s.mrcs' % (T, name)
    if not os.path.exists(fn): raise SystemExit('merge55: missing stack ' + fn)
    h = np.fromfile(fn, dtype=np.int32, count=256); nx, ny, nz, nsym = int(h[0]), int(h[1]), int(h[2]), int(h[23])
    mm.append(np.memmap(fn, dtype=np.float32, mode='r', offset=1024 + nsym, shape=(nz, ny, nx)))
    o = [x.rstrip('\n') for x in open('%s/oris_%s.txt' % (T, name)) if x.strip()]; assert len(o) == nz, (name, nz, len(o))
    oris += o; lab += [l] * nz; src.append(np.stack([np.full(nz, k), np.arange(nz)], 1))
lab = np.array(lab); src = np.concatenate(src); N = len(lab); perm = np.random.default_rng(11).permutation(N)
hdr = np.fromfile('%s/stk_%s.mrcs' % (T, names[0]), dtype=np.int32, count=256).copy(); hdr[2] = N; hdr[9] = N; hdr[23] = 0
CH = int(os.environ.get('CHUNK', '20000'))
with open(T + '/all.mrcs', 'wb') as f:
    f.write(hdr.tobytes())
    for c0 in range(0, N, CH):
        idx = perm[c0:c0 + CH]; buf = np.empty((len(idx), ny, nx), np.float32)
        for k in np.unique(src[idx, 0]):
            w = np.where(src[idx, 0] == k)[0]; loc = src[idx[w], 1]; o = np.argsort(loc); buf[w[o]] = mm[k][loc[o]]
        f.write(buf.tobytes())
with open(T + '/all.txt', 'w') as fa, open(T + '/imp.txt', 'w') as fi, open(T + '/gt_particles.txt', 'w') as fg:
    fg.write('# particle conformer closed e1 e2 (truth; conformer 0 = junk)\n')
    for j, i in enumerate(perm):
        d = kv(oris[i]); d['state'] = str(lab[i]) if lab[i] > 0 else '0'
        fa.write(' '.join('%s=%s' % (k, v) for k, v in d.items()) + '\n')
        fi.write(' '.join('%s=%s' % (k, d[k]) for k in ('dfx', 'dfy', 'angast', 'phshift') if k in d) + '\n')
        fg.write('%d %d %d %s %s\n' % (j + 1, lab[i], int(lab[i] in closed), d.get('e1', '0'), d.get('e2', '0')))
print('merged %d particles: %d conformer, %d junk; closed %d' % (N, (lab > 0).sum(), (lab == 0).sum(), sum(int(l in closed) for l in lab)))

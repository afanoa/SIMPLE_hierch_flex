#!/usr/bin/env python3
"""Split lattice of a hetero_tree: every node's stage-B state pairs give child-difference maps (even / odd unfiltered halves, low-passed,
inside an envelope from the root mean); splits whose difference maps agree across halves (|cross-half cosine| >= TAU, even of one x odd
of the other: no shared noise) are the same motion found in different branches and are tied into one class (single linkage) = a factor
of the lattice. No spatial basis: the classes are made of the tree's own splits. Truth (optional) only scores: each class's energy share
in the truth domain masks.
usage: mm_lattice.py TREE_ROOT(1_hetero_tree copy with node_*/stageB/recvol_stateNN_{even,odd}_unfil.mrc) [tau=0.9] [lp=12] [masks=DIR]"""
import sys, os, glob, itertools, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
from scipy.ndimage import zoom, gaussian_filter, binary_dilation
args = dict(a.split('=', 1) if '=' in a else ('tree', a) for a in sys.argv[1:]); T = args['tree']; tau = float(args.get('tau', 0.9)); LP = float(args.get('lp', 12))
nodes = sorted({p.split('/node_')[1].split('/')[0] for p in glob.glob(T + '/node_r*/stageB/recvol_state*_even_unfil.mrc')}, key=lambda x: (len(x), x))
def ld(nd, s, h): return read_mrc('%s/node_%s/stageB/recvol_state%02d_%s_unfil.mrc' % (T, nd, s, h)).astype(np.float64)
B = ld(nodes[0], 1, 'even').shape[0]; smpd = 160 * 1.464 / B
f = np.fft.fftfreq(B, smpd); fr = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2); filt = fr <= 1. / LP
lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * filt))
mean = lp(ld('r', 1, 'even') + ld('r', 1, 'odd')); sm = gaussian_filter(mean, 2); env = binary_dilation(sm > 0.15 * sm.max(), iterations=4)
M = []
if 'masks' in args:
    k = 1
    while os.path.exists('%s/domain%d_mask.mrc' % (args['masks'], k)):
        d = read_mrc('%s/domain%d_mask.mrc' % (args['masks'], k)).astype(np.float32); M.append((zoom(d, B / d.shape[0], order=1) if d.shape[0] != B else d) > 0.5); k += 1
dirs = []
for nd in nodes:
    ns = len(glob.glob('%s/node_%s/stageB/recvol_state*_even_unfil.mrc' % (T, nd)))
    if ns < 2: continue
    E = [lp(ld(nd, s, 'even'))[env] for s in range(1, ns + 1)]; O = [lp(ld(nd, s, 'odd'))[env] for s in range(1, ns + 1)]
    for a, b in itertools.combinations(range(ns), 2):
        de, do = E[a] - E[b], O[a] - O[b]; nn = de @ do
        if nn <= 0: continue   # a difference that does not reproduce across halves is noise
        e = np.zeros_like(mean); e[env] = (0.5 * (de + do)) ** 2
        dirs.append(('%s:%d-%d' % (nd, a + 1, b + 1), de, do, np.sqrt(nn), ' '.join('d%d %.2f' % (k + 1, e[m].sum() / e.sum()) for k, m in enumerate(M))))
n = len(dirs); Cs = np.eye(n)
for i, j in itertools.combinations(range(n), 2):
    Cs[i, j] = Cs[j, i] = 0.5 * (dirs[i][1] @ dirs[j][2] + dirs[j][1] @ dirs[i][2]) / (dirs[i][3] * dirs[j][3])
cls = -np.ones(n, int); c = 0
for i in range(n):
    if cls[i] >= 0: continue
    st = [i]; cls[i] = c
    while st:
        u = st.pop()
        for v in np.where(np.abs(Cs[u]) >= tau)[0]:
            if cls[v] < 0: cls[v] = c; st.append(v)
    c += 1
print('%d split differences over %d nodes; tau %.2f -> %d classes' % (n, len(nodes), tau, c))
for k in range(c):
    mem = np.where(cls == k)[0]
    if len(mem) < 2: continue
    print('class %d (%d splits, %d nodes): %s' % (k, len(mem), len({dirs[i][0].split(':')[0] for i in mem}), '; '.join('%s [%s]' % (dirs[i][0], dirs[i][4]) for i in mem)))
single = [dirs[i][0] + (' [%s]' % dirs[i][4] if dirs[i][4] else '') for i in range(n) if (cls == cls[i]).sum() == 1]
print('untied splits (%d): %s' % (len(single), '; '.join(single)))

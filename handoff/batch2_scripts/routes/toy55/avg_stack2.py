#!/usr/bin/env python3
"""Stitching input: the conditional sub-class averages of a cls_expansion run (full, even-half and odd-half models) as
one particle stack with one oritab row per average carrying the sub-class's member-weighted projection direction
(e1 e2 from the refined ptcl3D poses; e3 = x = y = 0 because the averages live in the class frame, which for
direction-bin classes is the 3D frame at the bin direction). Rows 1..N = full averages, N+1..2N even, 2N+1..3N odd.
usage: avg_stack.py CLSEXP_DIR ptcl3D_refined.txt OUT_DIR"""
import sys, os, numpy as np
d, oris, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
def read_stack(fn):
    h = np.frombuffer(open(fn, 'rb').read(1024), np.int32); nx, ny, nz, nsym = [int(v) for v in (h[0], h[1], h[2], h[23])]
    return np.array(np.memmap(fn, np.float32, 'r', offset=1024 + nsym, shape=(nz, ny, nx)), np.float32), bytearray(open(fn, 'rb').read(1024))
full, hdr = read_stack(os.path.join(d, 'cls_expansion_cavgs.mrc')); even, _ = read_stack(os.path.join(d, 'cls_expansion_cavgs_even.mrc')); odd, _ = read_stack(os.path.join(d, 'cls_expansion_cavgs_odd.mrc'))
N = full.shape[0]; assert even.shape[0] == N and odd.shape[0] == N, 'half stacks differ in length'
lab = np.loadtxt(os.path.join(d, 'cls_expansion_labels.txt'), comments='#').astype(int)   # particle_index parent_class global_subclass
cmap = np.loadtxt(os.path.join(d, 'cls_expansion_class_map.txt'), comments='#')  # global_subclass parent_class local_subclass pop ...
gids = cmap[:, 0].astype(int); assert len(gids) == N, 'class map rows (%d) != stack images (%d)' % (len(gids), N); row_of = {g: j for j, g in enumerate(gids)}   # stack row j holds global sub-class gids[j] (any order)
D = {}
for ln in open(oris):
    if 'e1=' not in ln: continue
    kv = dict(p.split('=', 1) for p in ln.split() if '=' in p); D.setdefault('e1', []).append(float(kv['e1'])); D['e2'] = D.get('e2', []) + [float(kv['e2'])]
e1 = np.array(D['e1']); e2 = np.array(D['e2']); v = np.c_[np.sin(np.radians(e2)) * np.cos(np.radians(e1)), np.sin(np.radians(e2)) * np.sin(np.radians(e1)), np.cos(np.radians(e2))]
dirs = np.zeros((N, 3)); pop = np.zeros(N, int); parent = cmap[:, 1].astype(int)
for j, g in enumerate(gids):
    m = lab[lab[:, 2] == g, 0] - 1; pop[j] = len(m)
    if len(m): u = v[m].mean(0); dirs[j] = u / np.linalg.norm(u)
    else: dirs[j] = [0, 0, 1]
spread = np.array([np.degrees(np.arccos(np.clip(v[lab[lab[:, 2] == g, 0] - 1] @ dirs[j], -1, 1))).mean() if pop[j] else 0 for j, g in enumerate(gids)])
E2 = np.degrees(np.arccos(np.clip(dirs[:, 2], -1, 1))); E1 = np.degrees(np.arctan2(dirs[:, 1], dirs[:, 0])) % 360.
print('%d sub-class averages, members %d-%d, mean angular spread about the sub-class direction %.1f deg (max %.1f)' % (N, pop.min(), pop.max(), spread.mean(), spread.max()))
# stack: per-image normalisation (zero mean, unit sd inside the box) to a common scale
S = np.concatenate([full, even, odd], 0).astype(np.float32); S = (S - S.mean((1, 2), keepdims=True)) / (S.std((1, 2), keepdims=True) + 1e-12)
ii = np.frombuffer(bytes(hdr), np.int32).copy(); ii[2] = 3 * N; ii[9] = 3 * N; ii[23] = 0; ff = ii.view(np.float32); ff[19] = float(S.min()); ff[20] = float(S.max()); ff[21] = float(S.mean())
open(os.path.join(out, 'avgs.mrcs'), 'wb').write(ii.tobytes() + S.tobytes())
with open(os.path.join(out, 'avgs_oritab.txt'), 'w') as f:
    for rep in range(3):
        for g in range(N): f.write('e1=%.4f e2=%.4f e3=0.0000 x=0.0000 y=0.0000 state=1\n' % (E1[g], E2[g]))
np.savez(os.path.join(out, 'avgs_meta.npz'), N=N, gid=gids, parent=parent, pop=pop, dirs=dirs, spread=spread, e1=E1, e2=E2)   # consumers map labels' global_subclass -> stack row via gid
print('wrote %s/avgs.mrcs (%d images) and avgs_oritab.txt' % (out, 3 * N))

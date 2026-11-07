#!/usr/bin/env python3
"""Lattice start positions for the Fortran board from a latent (the prototype's --init npz / pca12): the first two
standardised columns of the meta's latent, 5th-95th percentiles mapped onto [0, n1-1] x [0, n2-1], one "pind g1 g2"
line per particle of the subset. usage: klem_init_positions.py META SUBSET OUT [--n1 7 --n2 7 --zkey z692]"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('meta'); ap.add_argument('subset'); ap.add_argument('out'); ap.add_argument('--n1', type=int, default=7); ap.add_argument('--n2', type=int, default=7); ap.add_argument('--zkey', default='')
a = ap.parse_args(); M = np.load(a.meta); keys = [k for k in M.keys() if k.startswith('z')]; zk = a.zkey or (keys[0] if keys else None)
assert zk is not None, 'no latent in the meta'
Z = M[zk].astype(float); flags = M['flags'].astype(int); N = len(flags); ref = np.where(flags == 1)[0]
pos = -np.ones(N, int); pos[ref] = np.arange(len(ref)) if len(Z) == len(ref) else np.arange(N)[ref]
sub = np.loadtxt(a.subset, dtype=int); rows = sub - 1; ok = pos[rows] >= 0 if len(Z) == len(ref) else np.ones(len(rows), bool)
zs = (Z - Z.mean(0)) / np.maximum(Z.std(0), 1e-12); zz = zs[pos[rows[ok]]] if len(Z) == len(ref) else zs[rows[ok]]
def sc(v, nn): return np.clip((v - np.percentile(v, 5)) / max(np.percentile(v, 95) - np.percentile(v, 5), 1e-9), 0, 1) * (nn - 1)
g1 = sc(zz[:, 0], a.n1); g2 = sc(zz[:, 1], a.n2)
np.savetxt(a.out, np.c_[sub[ok], g1, g2], fmt='%d %.4f %.4f'); print('%s: %d of %d particles placed from %s PCs 1-2 (std %s)' % (a.out, ok.sum(), len(sub), zk, np.array2string(Z.std(0)[:4], precision=3)))

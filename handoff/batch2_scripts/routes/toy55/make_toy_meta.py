#!/usr/bin/env python3
"""A meta npz for a toy project in the real-data chain's format: flags (1 for every row), poses and CTF from a ptcl3D
table (print_project_field), eo by row parity, stkind/indstk (one stack, row order), and a stand-in latent z692 (random)
so the lattice init has something to read. usage: make_toy_meta.py PTCL3D_TABLE OUT.npz"""
import sys, numpy as np
tab, out = sys.argv[1:3]; rows = []
for l in open(tab):
    if not l.strip() or l.startswith('#'): continue
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); rows.append([float(kv.get(k, 0)) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')])
A = np.array(rows); N = len(A); rng = np.random.default_rng(0)
np.savez(out, flags=np.ones(N, np.int32), e1=A[:, 0], e2=A[:, 1], e3=A[:, 2], x=A[:, 3], y=A[:, 4], dfx=A[:, 5], dfy=A[:, 6], angast=A[:, 7], eo=(np.arange(N) % 2).astype(np.int32),
         stkind=np.ones(N, np.int32), indstk=np.arange(1, N + 1, dtype=np.int32), z692=rng.standard_normal((N, 4)).astype(np.float32))
print('%s: %d rows (flags all 1, random stand-in latent)' % (out, N))

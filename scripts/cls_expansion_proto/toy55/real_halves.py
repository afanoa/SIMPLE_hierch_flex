#!/usr/bin/env python3
"""Gold-standard halves of the real RNase P refined set for the data-derived tree: rows with flags == 1 split by the
refinement's even/odd field; per half a row-ordered stack (every image read through its own stkind / indstk keys, never
by row number: the project rows are a permutation of the stack order), the pose table (print_project_field style), a
placeholder truth file (conformer 1 for every particle: real data has no truth) and the meta row of every half row.
One pass per source stack, both halves filled through write memmaps. usage: real_halves.py META_NPZ STACK_GLOB SIM_PREFIX LAT_PREFIX"""
import sys, os, glob, numpy as np
meta, sglob, simp, latp = sys.argv[1:5]
M = {k: v for k, v in np.load(meta).items() if k != 'z692'}; f = M['flags'].astype(int); eo = M['eo'].astype(int); stk = M['stkind'].astype(int); ind = M['indstk'].astype(int)
stacks = sorted(glob.glob(sglob)); assert len(stacks) == stk.max(), (len(stacks), stk.max())
h0 = np.fromfile(stacks[0], np.int32, 256).copy(); nx, ny = int(h0[0]), int(h0[1])
halves = {}
for half, e in (('A', 0), ('B', 1)):
    rows = np.where((f == 1) & (eo == e))[0]; N = len(rows); d = '%s_%s' % (simp, half); os.makedirs(d, exist_ok=True); os.makedirs('%s_%s' % (latp, half), exist_ok=True)
    hdr = h0.copy(); hdr[2] = N; hdr[9] = N; hdr[23] = 0; open(d + '/all.mrcs', 'wb').write(hdr.tobytes())
    with open(d + '/all.mrcs', 'r+b') as fo: fo.truncate(1024 + N * nx * ny * 4)
    halves[half] = (rows, np.memmap(d + '/all.mrcs', np.float32, 'r+', offset=1024, shape=(N, ny, nx)))
    np.savetxt(d + '/meta_rows.txt', rows + 1, fmt='%d')
    with open('%s_%s/ptcl3D_r0.txt' % (latp, half), 'w') as fo:
        for i in rows: fo.write('e1=%.4f e2=%.4f e3=%.4f x=%.4f y=%.4f dfx=%.5f dfy=%.5f angast=%.3f state=1\n' % tuple(float(M[k][i]) for k in ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast')))
    with open(d + '/gt_particles.txt', 'w') as fo:
        fo.write('# particle conformer closed e1 e2 (REAL DATA: no truth, conformer 1 for every particle)\n')
        for j, i in enumerate(rows): fo.write('%d 1 0 %.3f %.3f\n' % (j + 1, M['e1'][i], M['e2'][i]))
    print('half %s: %d particles (eo %d)' % (half, N, e), flush=True)
for s, fn in enumerate(stacks, 1):
    h = np.fromfile(fn, np.int32, 256); nz, nsym = int(h[2]), int(h[23]); assert int(h[3]) == 2 and int(h[0]) == nx
    S = np.fromfile(fn, np.float32, offset=1024 + nsym).reshape(nz, ny, nx)
    for half, (rows, out) in halves.items():
        m = np.where(stk[rows] == s)[0]; out[m] = S[ind[rows[m]] - 1]
    print('stack %d (%s): %d images placed' % (s, os.path.basename(fn), sum(int((stk[r] == s).sum()) for r, _ in halves.values())), flush=True); del S
for half, (rows, out) in halves.items(): out.flush(); print('half %s stack mean %.4f sd %.4f (first 2000)' % (half, float(out[:2000].mean()), float(out[:2000].std())))

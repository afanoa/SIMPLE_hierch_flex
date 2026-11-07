#!/usr/bin/env python3
"""Data for the f692rf tree / r22 UMAP figures: (1) tree.json: f692rf nodes (queue.txt) and leaves (leaves.txt) with particles, the old
f692r closed states 15+16 traced in (count, share), the leaves' closed-arm / open-arm readout, plus the unit-7 extension (rflex_f692u7);
(2) r22_umap.npz: UMAP of r22's flex latent (centred per direction bin, standardised; NSUB subsample) with per particle the r22 refined
state, the unit-7 refined state (0 outside unit 7), old-closed flag and P(closed site occupied) from site_r22.
usage: viz_prep.py OUTDIR [--nsub 60000]"""
import os, re, json, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--nsub', type=int, default=60000); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'; T = '/mnt/beegfs/elmlund/afan/routes/toy55'; R = K_ + '/rflex_f692rf'; U7 = K_ + '/rflex_f692u7/node_r'
old = np.loadtxt(K_ + '/f692r/coordtree_fam1/units.txt', dtype=int); cl = (old == 15) | (old == 16); NR = len(old)
def arm_file(p):
    d = {}
    if os.path.isfile(p):
        for l in open(p):
            if l.startswith('#') or not l.strip(): continue
            f = l.split(); d[os.path.basename(f[0]).replace('.mrc', '')] = (float(f[1]), float(f[2]))
    return d
nodes = []
for l in open(R + '/queue.txt'):
    f = l.split()
    if len(f) < 4: continue
    rows = np.loadtxt(f[3]) > 0; nodes.append(dict(id=f[0], depth=int(f[1]), n=int(rows.sum()), closed=int(cl[rows].sum())))
arm = arm_file(K_ + '/view_rflex/info/arm_f692rf_leaves.txt')
for i, l in enumerate(open(R + '/leaves.txt'), 1):
    f = l.split()
    if len(f) < 7: continue
    lab = np.loadtxt(f[5], dtype=int); m = lab == int(f[6]); nm = f[0]; parent = nm.rsplit('_', 1)[0]
    key = [k for k in arm if k.startswith('leaf%02d_' % i)]; ca = arm[key[0]] if key else (np.nan, np.nan)
    nodes.append(dict(id=nm, parent=parent, leaf=True, n=int(m.sum()), closed=int(cl[m].sum()), closed_arm=ca[0], open_arm=ca[1], junk=('junk' in l)))
lab7 = np.loadtxt(U7 + '/refined_f.txt', dtype=int); a7 = arm_file(U7 + '/readout_fstageB/arm.txt')
for s in range(1, lab7.max() + 1):
    m = lab7 == s; ca = a7.get('state_%02d' % s, (np.nan, np.nan))
    nodes.append(dict(id='u7_s%d' % s, parent='r22_s2', leaf=True, ext=True, n=int(m.sum()), closed=int(cl[m].sum()), closed_arm=ca[0], open_arm=ca[1]))
json.dump(dict(nodes=nodes, total_closed=int(cl.sum())), open(a.out + '/tree.json', 'w'), indent=1); print('tree.json: %d entries' % len(nodes))
# ---- r22 UMAP
F = np.loadtxt(R + '/node_r22/1_flex_pca/flex_pca_coordinates.txt', comments='#', ndmin=2); pidx = F[:, 0].astype(int); Z = F[:, 5:]
e1 = np.zeros(NR); e2 = np.zeros(NR); i = 0
for ln in open(K_ + '/cheapkids_f692rf_r22/ptcl3D_parent.txt'):
    if 'e1=' not in ln: continue
    kv = dict(t.split('=', 1) for t in ln.split() if '=' in t); e1[i] = float(kv['e1']); e2[i] = float(kv['e2']); i += 1
rng = np.random.default_rng(0); E1, E2 = np.radians(e1[pidx - 1]), np.radians(e2[pidx - 1]); v = np.c_[np.sin(E2) * np.cos(E1), np.sin(E2) * np.sin(E1), np.cos(E2)]
N = len(pidx); nb = max(20, N // 1300); C = v[rng.choice(N, nb, replace=False)].copy()
for it in range(30):
    lb = np.argmax(v @ C.T, 1)
    for k in range(nb):
        mm = lb == k; C[k] = v[mm].mean(0) / np.linalg.norm(v[mm].mean(0)) if mm.any() else v[rng.integers(N)]
lb = np.argmax(v @ C.T, 1); _, bi = np.unique(lb, return_inverse=True); cnt = np.bincount(bi)
bm = np.stack([np.bincount(bi, Z[:, j], minlength=bi.max() + 1) for j in range(Z.shape[1])], 1) / cnt[:, None]; Zc = Z - bm[bi]; Zc = (Zc - Zc.mean(0)) / Zc.std(0)
sub = np.sort(rng.choice(N, min(a.nsub, N), replace=False)); rows = pidx[sub] - 1
import umap
Y = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0).fit_transform(Zc[sub])
st22 = np.loadtxt(R + '/node_r22/refined_states.txt', dtype=int); st7 = np.zeros(NR, int); st7[:len(lab7)] = lab7
pc = np.full(NR, np.nan); H = T + '/simrnasep692'
for h in 'AB':
    mr = np.loadtxt(H + '_%s/meta_rows.txt' % h).astype(int) - 1; P = np.loadtxt(K_ + '/site_r22/em/post_soft_%s.txt' % h)[:, 1]
    m = np.loadtxt(K_ + '/site_r22/start_%s.txt' % h, dtype=int) > 0; pc[mr[m]] = P[m]
np.savez_compressed(a.out + '/r22_umap.npz', xy=Y.astype(np.float32), state=st22[rows], u7=st7[rows], oldclosed=cl[rows], pclosed=pc[rows].astype(np.float32), row=rows)
print('r22 UMAP: %d of %d particles, %d direction bins; old closed in the subsample %d; P(closed site) available for %d' % (len(sub), N, bi.max() + 1, cl[rows].sum(), np.isfinite(pc[rows]).sum()))

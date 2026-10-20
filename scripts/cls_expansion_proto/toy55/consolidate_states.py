#!/usr/bin/env python3
"""Consolidated states of a run: cut the leaf merge tree (merge_states.py) at a similarity, pool the particles of every
merged group (union over the member leaf states' selected lists), write OUT/states.txt (one merged-state label per project
row, 0 = none; a particle in several groups goes to the group where it has the most members) and OUT/groups.txt.
usage: consolidate_states.py RUNDIR MERGE.npz CUT NROWS OUT"""
import sys, os, re, numpy as np
from scipy.cluster.hierarchy import fcluster
run, mz, cut, nrows, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; os.makedirs(out, exist_ok=True)
M = np.load(mz, allow_pickle=True); names, Z = list(M['names']), (M['Z'] if 'Z' in M.files else None); lab = M['labels'].astype(int) if 'labels' in M.files else fcluster(Z, t=1 - cut, criterion='distance'); K = lab.max()
votes = np.zeros((nrows, K), np.float32); uncap = os.environ.get('UNCAP', '0') == '1'; _lab = {}
def labels_of(f):   # UNCAP=1: every particle whose argmax cell belongs to the state (clusters.txt cells x board latent), not the capped top-N selection
    if not uncap: st = np.loadtxt(f, dtype=int); return st, np.ones(len(st), np.float32)
    L1 = os.path.dirname(f)
    if L1 not in _lab:
        lat = np.loadtxt(L1 + '/../board/klem_latent.txt', comments='#'); pind, am = lat[:, 0].astype(int), lat[:, 3].astype(int); cell2state = {}
        for l in open(L1 + '/clusters.txt'):
            if l.startswith('#') or not l.strip(): continue
            t = l.split(); s_ = int(t[0]); [cell2state.__setitem__(int(c), s_) for c in t[4:]]
        lab = np.zeros(nrows, int); lab[pind - 1] = np.array([cell2state.get(c, 0) for c in am]); w = np.zeros(nrows, np.float32); w[pind - 1] = lat[:, 4]; _lab[L1] = (lab, w)
    return _lab[L1]
for i, nm in enumerate(names):
    if lab[i] == 0: continue   # label 0 = dropped leaf (cross-tree merge: unpaired)
    m = re.match(r'f(\d)_(.*)_s(\d+)$', nm); fam, reg, s = m.group(1), m.group(2), int(m.group(3))
    f = '%s/fam%s/region%s/level1/states_clusters.txt' % (run, fam, reg); f = f if os.path.isfile(f) else f.replace('/region', '/zoomregion'); st, w = labels_of(f); rows = np.where(st == s)[0]; votes[rows, lab[i] - 1] += w[rows]   # every region board labels every particle: ties are broken by the board's confidence (maxr), not by group index
states = np.where(votes.sum(1) > 0, votes.argmax(1) + 1, 0); sizes = np.bincount(states, minlength=K + 1)[1:]
minp = int(os.environ.get('MINP', '200')); small = np.where(sizes < minp)[0] + 1   # groups below MINP particles are not states (and crash the reconstruction)
if len(small): states[np.isin(states, small)] = 0; print('dropped %d merged groups below %d particles: %s' % (len(small), minp, ' '.join(str(sizes[c - 1]) for c in small)))
keepc = np.where(sizes >= minp)[0]; remap0 = np.zeros(K + 1, int); remap0[keepc + 1] = np.arange(1, len(keepc) + 1); states = remap0[states]; lab = np.array([remap0[l] for l in lab]); K = len(keepc); sizes = np.bincount(states, minlength=K + 1)[1:]
order = np.argsort(-sizes); remap = np.zeros(K + 1, int); remap[order + 1] = np.arange(1, K + 1); states = remap[states]
if os.environ.get('VIEWMATCH', '0') == '1':   # view-matched states: thin every state to the whole set's direction histogram (100 bins), so that two states of one conformer are not told apart by view coverage
    meta = os.environ.get('META', run + '/fam1/meta_fb.npz'); Mt = np.load(meta); e1, e2 = np.radians(Mt['e1']), np.radians(Mt['e2']); v = np.c_[np.sin(e2) * np.cos(e1), np.sin(e2) * np.sin(e1), np.cos(e2)]
    rng = np.random.default_rng(0); C = v[rng.choice(len(v), 100, replace=False)].copy()
    for it in range(30):
        vb = np.argmax(v @ C.T, 1)
        for k in range(100):
            mk = vb == k
            if mk.any(): C[k] = v[mk].mean(0); C[k] /= np.linalg.norm(C[k])
    vb = np.argmax(v @ C.T, 1)[:nrows]; H0 = np.bincount(vb[states > 0], minlength=100) / max((states > 0).sum(), 1); kept = 0; ntot = int((states > 0).sum())
    for s in range(1, K + 1):
        rows_s = np.where(states == s)[0]
        if len(rows_s) == 0: continue
        H = np.bincount(vb[rows_s], minlength=100) / len(rows_s); keep_p = np.minimum(1., H0[vb[rows_s]] / np.maximum(H[vb[rows_s]], 1e-9)); drop = rng.random(len(rows_s)) > keep_p; states[rows_s[drop]] = 0; kept += int((~drop).sum())
    print('view matching: %d particles kept of %d (states thinned to the whole set\'s direction histogram)' % (kept, ntot))
np.savetxt(out + '/states.txt', states, fmt='%d')
with open(out + '/groups.txt', 'w') as o:
    o.write('# merged state (cut %.2f): particles, members (leaf states)\n' % cut)
    for r, c in enumerate(order): o.write('%d %d %s\n' % (r + 1, sizes[c], ' '.join(names[i] for i in range(len(names)) if lab[i] == c + 1)))
print('cut %.2f: %d merged states from %d leaf states (%s); particles per merged state: %s' % (cut, K, len(names), 'ALL particles of each state' if uncap else 'capped leaf selections', ' '.join(str(sizes[c]) for c in order)))

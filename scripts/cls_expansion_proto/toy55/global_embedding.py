#!/usr/bin/env python3
"""ONE embedding of the whole tree: every particle is represented by the concatenation of its responsibility vectors
over all the boards it took part in (parent board, its family's feedback board, every region board of the family, the
level-2 boards it fell into), Hellinger-transformed (sqrt), and UMAP-embedded jointly, so the map reflects every view
the workflow produced at once. Islands are counted by HDBSCAN (fallback DBSCAN) on the embedding, with the dominant
family and identified conformer per island. Colourings: family, identified conformer (feedback-board cell match),
junk mass, defocus, and the count_states landscape state when a leaf clustering file is given.
usage: global_embedding.py RUNDIR META OUT.png [--nsub 80000 --min_island 500 --seed 0]"""
import sys, os, glob, argparse, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('meta'); ap.add_argument('out'); ap.add_argument('--nsub', type=int, default=80000); ap.add_argument('--min_island', type=int, default=500); ap.add_argument('--seed', type=int, default=0); ap.add_argument('--gt', default='', help='toy truth: gt_particles.txt (particle conformer closed ...)'); ap.add_argument('--holo0', type=int, default=101)
a = ap.parse_args(); R = a.run; M = np.load(a.meta); N = len(M['flags']); rng = np.random.default_rng(a.seed)
def load_resp(d):
    if not os.path.exists(d + '/klem_resp.npz'): return None
    Z = np.load(d + '/klem_resp.npz'); return Z['pinds'].astype(int), Z['r'].astype(np.float32)
boards = [('parent', R + '/parent')]
for fam in sorted(glob.glob(R + '/fam?')):
    f = os.path.basename(fam); boards.append((f + '/feedback', fam + '/klem_fb'))
    for rb in sorted(glob.glob(fam + '/region*/board')): boards.append((f + '/' + rb.split('/')[-2], rb))
    for l2 in sorted(glob.glob(fam + '/region*/level2_s*/board')): boards.append((f + '/' + l2.split('/')[-3] + '/' + l2.split('/')[-2], l2))
loaded = [(name, load_resp(d)) for name, d in boards]; loaded = [(n, v) for n, v in loaded if v is not None]
print('boards in the embedding: %d' % len(loaded)); [print('  %-40s %7d particles x %d components' % (n, len(v[0]), v[1].shape[1])) for n, v in loaded]
# particle universe = the parent board's on-grid particles; subsample
pp, rp = loaded[0][1]; universe = pp; sub = rng.choice(len(universe), min(a.nsub, len(universe)), replace=False); sub.sort(); pinds = universe[sub]
pos = {p: i for i, p in enumerate(pinds)}; blocks = []; cols = 0
for name, (pi, r) in loaded:
    blk = np.zeros((len(pinds), r.shape[1]), np.float32); idx = np.array([pos.get(p, -1) for p in pi]); m = idx >= 0; blk[idx[m]] = r[m]
    blocks.append(np.sqrt(np.maximum(blk, 0))); cols += r.shape[1]
X = np.concatenate(blocks, 1); print('feature matrix: %d particles x %d (sqrt responsibilities over all boards)' % X.shape)
import umap
emb = umap.UMAP(n_neighbors=30, min_dist=0.05, metric='euclidean', random_state=a.seed).fit_transform(X)
# islands
try:
    import hdbscan; lab = hdbscan.HDBSCAN(min_cluster_size=a.min_island).fit_predict(emb); how = 'HDBSCAN'
except Exception:
    from sklearn.cluster import DBSCAN; span = float(np.max(emb.max(0) - emb.min(0))); lab = DBSCAN(eps=0.02 * span, min_samples=30).fit_predict(emb); how = 'DBSCAN(eps 2%% of the span)'
isl = [k for k in np.unique(lab) if k >= 0 and (lab == k).sum() >= a.min_island]
# annotations: family, identified conformer (feedback-board cell match), junk mass (parent free components), defocus
fam_of = np.zeros(len(pinds), int); conf_of = np.zeros(len(pinds), int)
for j, fam in enumerate(sorted(glob.glob(R + '/fam?'))):
    s = np.loadtxt(R + '/expand/family%d_subset.txt' % (j + 1), dtype=int); fs = set(s.tolist()); fam_of[[i for i, p in enumerate(pinds) if p in fs]] = j + 1
    cm = fam + '/klem_fb/cell_matches.txt'
    if os.path.exists(cm) and os.path.exists(fam + '/klem_fb/klem_latent.txt'):
        Mt = np.loadtxt(cm, comments='#'); best = {int(r[0]): int(r[1]) for r in Mt}; L = np.loadtxt(fam + '/klem_fb/klem_latent.txt', comments='#')
        for p, c, g in zip(L[:, 0].astype(int), L[:, 3].astype(int), L[:, 5]):
            if p in pos and g > 0.5 and c in best: conf_of[pos[p]] = best[c]
ngrid_p = len(np.load(R + '/parent/klem_resp.npz')['i1']); junk = 1 - rp[sub, :ngrid_p].sum(1); dfx = M['dfx'][pinds - 1] if 'dfx' in M else np.zeros(len(pinds))
print('%s islands (>= %d particles): %d; particles in islands %.2f, noise %.2f' % (how, a.min_island, len(isl), np.isin(lab, isl).mean(), (lab < 0).mean()))
for k in isl:
    m = lab == k; fams = np.bincount(fam_of[m]); cf = conf_of[m]; cf = cf[cf > 0]; top = np.bincount(cf).argmax() if len(cf) else 0
    print('  island %2d: %6d particles | family %s | identified conformer %s (%.0f%% of its labelled particles) | junk mass %.2f' % (k, m.sum(), fams.argmax(), top, 100 * (cf == top).mean() if len(cf) else 0, junk[m].mean()))
gt = None
if a.gt:
    g_ = {}
    for l in open(a.gt):
        if l.strip() and not l.startswith('#'): w = l.split(); g_[int(w[0])] = (int(w[1]), int(w[2]))
    tconf = np.array([g_.get(p, (0, 0))[0] for p in pinds]); tclosed = np.array([g_.get(p, (0, 0))[1] for p in pinds]); gt = (tconf, tclosed)
    print('TRUTH per island:')
    for k in isl:
        m = lab == k; c_ = tconf[m]; c_ = c_[c_ > 0]; top = np.bincount(c_).argmax() if len(c_) else 0
        print('  island %2d: %6d particles | true conformer %d (%.0f%% purity) | holo %.2f | closed %.2f | junk (true) %.2f' % (k, m.sum(), top, 100 * (c_ == top).mean() if len(c_) else 0, (c_ >= a.holo0).mean() if len(c_) else 0, tclosed[m][tconf[m] > 0].mean() if (tconf[m] > 0).any() else 0, (tconf[m] == 0).mean()))
    conformers = np.unique(tconf[tconf > 0]); print('  conformers with a majority island: %d of %d' % (sum(1 for c in conformers if any(np.bincount(tconf[lab == k][tconf[lab == k] > 0]).argmax() == c for k in isl if (tconf[lab == k] > 0).any())), len(conformers)))
fig, ax = plt.subplots(3 if gt is not None else 2, 3, figsize=(18, 16 if gt is not None else 11))
def sc(a_, c, title, cmap='tab20', s=1.5, **kw):
    o = np.argsort(np.random.default_rng(1).random(len(emb))); im = a_.scatter(emb[o, 0], emb[o, 1], c=np.asarray(c)[o], s=s, cmap=cmap, linewidths=0, **kw); a_.set_title(title, fontsize=10); a_.set_xticks([]); a_.set_yticks([]); return im
sc(ax[0, 0], fam_of, 'family (parent split)', cmap='Set1')
im = sc(ax[0, 1], np.where(conf_of > 0, conf_of, np.nan), 'deposited conformer matched by the feedback-board cell (readout, not truth)', cmap='tab20');
sc(ax[0, 2], lab, '%s islands: %d of >= %d particles (grey = unassigned)' % (how, len(isl), a.min_island), cmap='tab20')
sc(ax[1, 0], junk, 'junk mass (parent free components)', cmap='magma', vmin=0, vmax=1)
sc(ax[1, 1], dfx, 'defocus (nuisance check)', cmap='viridis')
ax[1, 2].hexbin(emb[:, 0], emb[:, 1], gridsize=90, cmap='Greys', bins='log'); ax[1, 2].set_title('density', fontsize=10); ax[1, 2].set_xticks([]); ax[1, 2].set_yticks([])
if gt is not None:
    sc(ax[2, 0], np.where(gt[0] > 0, gt[0], np.nan), 'TRUE conformer (toy)', cmap='tab20'); sc(ax[2, 1], np.where(gt[0] > 0, (gt[0] >= a.holo0).astype(float), np.nan), 'TRUE holo (1) / apo (0)', cmap='coolwarm', vmin=0, vmax=1); sc(ax[2, 2], np.where(gt[0] > 0, gt[1].astype(float), np.nan), 'TRUE closed (1) / open (0)', cmap='viridis', vmin=0, vmax=1)
fig.suptitle('%s: one embedding of the whole tree (%d particles, %d boards, %d responsibility components)' % (os.path.basename(R.rstrip('/')), len(pinds), len(loaded), cols), fontsize=11)
plt.tight_layout(); plt.savefig(a.out, dpi=120); np.savez_compressed(a.out.replace('.png', '.npz'), pinds=pinds, emb=emb, island=lab, family=fam_of, conformer=conf_of, junk=junk); print('->', a.out)

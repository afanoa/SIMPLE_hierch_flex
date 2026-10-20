#!/usr/bin/env python3
"""Gallery of every leaf state map of a klem run in the same view: per map one row = three orthogonal projections of the
10 A low-passed map (shared grey scale) and the projection of its deviation from the consensus (diverging scale, shared),
labelled with family/region/level/state, particle mass, reconstructed count and the best deposited match. Pages of NROW.
usage: volume_gallery.py RUNDIR CONSENSUS OUTPREFIX [--nrow 12] [--lp 10] [--smpd 1.464]"""
import sys, os, glob, re, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('cons'); ap.add_argument('out'); ap.add_argument('--nrow', type=int, default=12)
ap.add_argument('--lp', type=float, default=10.); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--sort', default='mass'); a = ap.parse_args()
def lowpass(v, lp):
    n = v.shape[0]; f = np.fft.fftfreq(n, d=a.smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2
    return np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-0.5 * k2 * (lp / 1.2)**2)))
def fsc_res(fe, fo, thr=0.143):   # FSC between the even/odd maps under a soft spherical mask, resolution at thr
    e, o = read_mrc(fe).astype(np.float32), read_mrc(fo).astype(np.float32); n = e.shape[0]; x = np.arange(n) - n // 2
    rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
    E, O = np.fft.fftn(e * msk), np.fft.fftn(o * msk); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k = np.sqrt(kx**2 + ky**2 + kz**2) * n
    sh = np.minimum(np.round(k).astype(int), n // 2); num = np.bincount(sh.ravel(), np.real(np.conj(E) * O).ravel()); d1 = np.bincount(sh.ravel(), np.abs(E).ravel() ** 2); d2 = np.bincount(sh.ravel(), np.abs(O).ravel() ** 2)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < thr)[0]; sres = (below[0] + 1) if len(below) else n // 2 - 1
    return n * a.smpd / max(sres, 1)
cons = lowpass(read_mrc(a.cons).astype(np.float32), a.lp)
rows = []
for d in sorted(glob.glob(a.run + '/fam?/' + os.environ.get('LEAFGLOB', 'region*'))):
    m = re.search(r'fam(\d)/(?:zoom)?region(.*)$', d); fam, reg = m.group(1), m.group(2)
    cl = {}; 
    if os.path.isfile(d + '/level1/clusters.txt'):
        for l in open(d + '/level1/clusters.txt'):
            if l.startswith('#'): continue
            p = l.split(); cl[int(p[0])] = (float(p[2]), int(p[3]))
    mt = {}
    if os.path.isfile(d + '/level1/state_matches.txt'):
        for l in open(d + '/level1/state_matches.txt'):
            if l.startswith('#'): continue
            p = l.split(); mt[int(p[0])] = (p[1], float(p[2]))
    cells = {}
    if os.path.isfile(d + '/level1/clusters.txt'):
        for l in open(d + '/level1/clusters.txt'):
            if l.startswith('#'): continue
            p = l.split(); cells[int(p[0])] = np.array([int(x) - 1 for x in p[4:]])
    H = np.load(d + '/board/klem_halves.npz') if os.path.isfile(d + '/board/klem_halves.npz') else None
    def xh(ca, cb):
        a0, a1, b0, b1 = H['D0'][ca].mean(0), H['D1'][ca].mean(0), H['D0'][cb].mean(0), H['D1'][cb].mean(0)
        cc = lambda x, y: np.real(np.sum(np.conj(x) * y)) / max(np.sqrt(np.sum(np.abs(x) ** 2) * np.sum(np.abs(y) ** 2)), 1e-30)
        return 0.5 * (cc(a0, b1) + cc(a1, b0))
    for f in sorted(glob.glob(d + '/level1/states_dir/state_*.mrc')):
        s = int(re.findall(r'state_(\d+)', f)[0]); mass, nsel = cl.get(s, (np.nan, 0)); c, cc = mt.get(s, ('?', np.nan))
        rep, sep = np.nan, np.nan
        rep2 = np.nan
        if H is not None and s in cells:
            rep = xh(cells[s], cells[s]); others = [xh(cells[s], cells[t]) for t in cells if t != s]; sep = max(others) if others else np.nan
            allc = np.concatenate(list(cells.values())); m0, m1 = H['D0'][allc].mean(0), H['D1'][allc].mean(0)
            e0, e1 = H['D0'][cells[s]].mean(0) - m0, H['D1'][cells[s]].mean(0) - m1
            rep2 = np.real(np.sum(np.conj(e0) * e1)) / max(np.sqrt(np.sum(np.abs(e0) ** 2) * np.sum(np.abs(e1) ** 2)), 1e-30)
        fe, fo = d + '/level1/maps/recvol_state%02d_even.mrc' % s, d + '/level1/maps/recvol_state%02d_odd.mrc' % s
        res = fsc_res(fe, fo) if os.path.isfile(fe) and os.path.isfile(fo) else np.nan
        rows.append(dict(f=f, label='f%s %s L1 s%02d' % (fam, reg, s), mass=mass, nsel=nsel, conf=c, corr=cc, rep=rep, sep=sep, res=res, rep2=rep2))
    for L in sorted(glob.glob(d + '/level2_s*')):
        s2 = int(re.findall(r'level2_s(\d+)', L)[0]); n2 = sum(1 for _ in open(L + '/states_top.txt')) if os.path.isfile(L + '/states_top.txt') else 0
        ntop = int(np.sum(np.loadtxt(L + '/states_top.txt') > 0)) if n2 else 0
        for f in sorted(glob.glob(L + '/states_dir/state_*.mrc')):
            fe, fo = L + '/maps/recvol_state01_even.mrc', L + '/maps/recvol_state01_odd.mrc'; res = fsc_res(fe, fo) if os.path.isfile(fe) and os.path.isfile(fo) else np.nan
            rows.append(dict(f=f, label='f%s %s L2 of s%02d' % (fam, reg, s2), mass=float(ntop), nsel=ntop, conf='-', corr=np.nan, rep=np.nan, sep=np.nan, res=res, rep2=np.nan))
if a.sort == 'mass': rows.sort(key=lambda r: -r['mass'] if np.isfinite(r['mass']) else 0)
elif a.sort == 'rep': rows.sort(key=lambda r: r['rep2'] if np.isfinite(r['rep2']) else 9)
print('%d maps' % len(rows))
vols = [lowpass(read_mrc(r['f']).astype(np.float32), a.lp) for r in rows]
supp = cons > 0.2 * cons.max(); outside = cons < 0.02 * cons.max()
def standardise(v):   # zero the solvent (mean outside the consensus envelope), unit RMS inside the support: every map on one scale
    v = v - v[outside].mean(); return v / max(np.sqrt(np.mean(v[supp] ** 2)), 1e-12)
cons = standardise(cons); vols = [standardise(v) for v in vols]
proj = [[v.sum(0), v.sum(1), v.sum(2)] for v in vols]; dev = [[(v - cons).sum(0), (v - cons).sum(1), (v - cons).sum(2)] for v in vols]
allp = np.concatenate([p.ravel() for pp in proj for p in pp]); vmin, vmax = np.percentile(allp, 1), np.percentile(allp, 99.8)
alld = np.concatenate([p.ravel() for pp in dev for p in pp]); dmax = np.percentile(np.abs(alld), 99.5)
npage = (len(rows) + a.nrow - 1) // a.nrow
for pg in range(npage):
    sel = range(pg * a.nrow, min(len(rows), (pg + 1) * a.nrow)); nr = len(sel)
    fig, axs = plt.subplots(nr, 6, figsize=(13.5, 2.1 * nr + 0.6), squeeze=False)
    for ri, i in enumerate(sel):
        r = rows[i]
        for j in range(3):
            axs[ri, j].imshow(proj[i][j].T, cmap='gray', vmin=vmin, vmax=vmax, origin='lower'); axs[ri, 3 + j].imshow(dev[i][j].T, cmap='RdBu_r', vmin=-dmax, vmax=dmax, origin='lower')
        for ax in axs[ri]: ax.set_xticks([]); ax.set_yticks([])
        junk = (np.isfinite(r['rep2']) and r['rep2'] < 0.5) or (np.isfinite(r['res']) and r['res'] > 20.)
        txt = '%s%s\nmass %s  rec %d  FSC %.1f A\nrep %.2f own %.2f sep %.2f\nbest %s r=%.2f' % ('WEAK? ' if junk else '', r['label'], ('%.0f' % r['mass']) if np.isfinite(r['mass']) else '?', r['nsel'], r['res'], r['rep'], r['rep2'], r['sep'], r['conf'], r['corr'])
        axs[ri, 0].set_ylabel(txt, rotation=0, ha='right', va='center', fontsize=8, labelpad=4)
    for j, t in enumerate(['proj z', 'proj y', 'proj x', 'dev z', 'dev y', 'dev x']): axs[0, j].set_title(t, fontsize=9)
    fig.suptitle('%s leaf state maps, %.0f A low-pass, page %d/%d (sorted by %s; dev = scale-matched map - consensus, shared scales)' % (os.path.basename(a.run.rstrip('/')), a.lp, pg + 1, npage, a.sort), fontsize=10)
    plt.subplots_adjust(left=0.22, right=0.99, top=0.95, bottom=0.01, wspace=0.03, hspace=0.08); fig.savefig('%s_p%02d.png' % (a.out, pg + 1), dpi=110); plt.close(fig)
with open(a.out + '_table.txt', 'w') as o:
    o.write('# rank label mass nrec fsc_A rep own_rep sep best_conformer corr file   (own_rep = cross-half reproducibility of the state deviation minus the region mean deviation)\n')
    for k, r in enumerate(rows): o.write('%d %s %s %d %.1f %.2f %.2f %.2f %s %.3f %s\n' % (k + 1, r['label'].replace(' ', '_'), ('%.0f' % r['mass']) if np.isfinite(r['mass']) else 'nan', r['nsel'], r['res'], r['rep'], r['rep2'], r['sep'], r['conf'], r['corr'], r['f']))
print('pages: %d -> %s_pNN.png, table %s_table.txt' % (npage, a.out, a.out))

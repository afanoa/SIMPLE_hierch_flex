#!/usr/bin/env python3
"""land_leafext.py: do the EXTREMES of independent runs' landscapes reproduce? Per run: its stage-B leaf half-maps
(MAPDIR/recvol_stateNN_{even,odd}_unfil.mrc) and leaf labels (one state per master row). As landscape_cmp.py: low-pass LP, node mask,
deviation maps from each run's leaf-size-weighted mean, half-map reliabilities of the deviations; a pooled leaf-size-weighted PCA over
all runs' leaves gives the stable directions. Per run and axis, its two most extreme leaves (lowest / highest projection). Across runs,
map against map, per pole:
  rho   = corr(d_A, d_B) / sqrt(R_A R_B)       normalised correlation of the deviation maps (1 = same within noise)
  delta = |F_A - F_B|^2 / (nu_A + nu_B)       noise-normalised difference of the full maps, nu = |E - O|^2 / 4 from each map's halves
                                              (1 = within noise for independent particles; shared particles push it below 1)
plus each pole's projection (pooled sd units) and the particle Jaccard of the two leaves. Controls: the same run's two poles, and
each pole against the other run's nearest non-pole leaf. --spec writes 'name even odd' lines of every extreme leaf for make_view.py.
usage: land_leafext.py NAME=MAPDIR:LABELS ... [--lp 15 --naxes 2 --mask M --spec OUT]"""
import sys, itertools, numpy as np
args = [x for x in sys.argv[1:] if '=' in x and not x.startswith('--')]; opt = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
LP = float(opt('--lp', 15.)); NAX = int(opt('--naxes', 2)); SPEC = opt('--spec', ''); FIG = opt('--fig', ''); K_ = '/mnt/beegfs/elmlund/afan/rnasep_klem'
MASK = opt('--mask', K_ + '/rflex_f692st/node_r1/rsmask.mrc')
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def cell(p): return float(np.fromfile(p, np.float32, 13)[10])
def crop(v, n):
    N = v.shape[0]
    if N == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); lo = N // 2 - n // 2; G = F[lo:lo + n, lo:lo + n, lo:lo + n]
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / N) ** 3
def lowpass(v, smpd, lp):
    n = v.shape[0]; f = np.fft.fftfreq(n, d=smpd); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij')
    return np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-0.5 * (kx ** 2 + ky ** 2 + kz ** 2) * (lp / 1.2) ** 2)))
RUNS = {}; n0 = None
for x in args:
    nm, rest = x.split('=', 1); mdir, labf = rest.split(':'); lab = np.loadtxt(labf, dtype=int); ks = sorted(set(lab[lab > 0].tolist())); E, O, files = [], [], []
    for k in ks:
        b = '%s/recvol_state%02d' % (mdir, k); e, o = rd(b + '_even_unfil.mrc'), rd(b + '_odd_unfil.mrc'); n = e.shape[0]; sm = cell(b + '_even_unfil.mrc') / n
        if n0 is None: n0, smpd0, cell0 = n, sm, cell(b + '_even_unfil.mrc')
        if n != n0: e, o = crop(e, n0), crop(o, n0)
        E.append(lowpass(e, smpd0, LP).ravel()); O.append(lowpass(o, smpd0, LP).ravel()); files.append((b + '_even_unfil.mrc', b + '_odd_unfil.mrc'))
    RUNS[nm] = dict(ks=ks, lab=lab, p=np.array([(lab == k).sum() for k in ks], float), E=np.array(E), O=np.array(O), files=files)
    print('%-7s %d leaves from %s (box %d)' % (nm, len(ks), mdir, n), flush=True)
Mk = np.clip(crop(rd(MASK), n0), 0, 1).ravel(); w = Mk / Mk.sum()
for R in RUNS.values():
    p = R['p'] / R['p'].sum(); F = 0.5 * (R['E'] + R['O']); dF, dE, dO = (X - (p[:, None] * X).sum(0) for X in (F, R['E'], R['O']))
    wc = lambda a, b: (w * a * b).sum() / np.sqrt((w * a * a).sum() * (w * b * b).sum())
    r = np.array([wc(dE[k], dO[k]) for k in range(len(p))]); R.update(pn=p, F=F, d=dF, rel=2 * r / (1 + r), nu=(w[None] * (R['E'] - R['O']) ** 2).sum(1) / 4.)
nr = len(RUNS); allD = np.vstack([R['d'] * np.sqrt(w)[None] * np.sqrt(R['pn'] / nr)[:, None] for R in RUNS.values()])
U, s, Vt = np.linalg.svd(allD, full_matrices=False); sh = s ** 2 / (s ** 2).sum(); V = Vt[:NAX]
print('mask %s, lp %.0f A: pooled leaf PCA variance shares %s' % (MASK.split('/')[-1], LP, ' '.join('%.3f' % x for x in sh[:4])))
for nm, R in RUNS.items():
    Dr = R['d'] * np.sqrt(w)[None] * np.sqrt(R['pn'])[:, None]; _, sr, Vr = np.linalg.svd(Dr, full_matrices=False)
    print('  %-7s own PC1/PC2 vs pooled: |cos| %s ; half-map reliability of its deviations %s' % (nm, ' '.join('%.2f' % abs(Vr[j] @ V[j]) for j in range(NAX)),
          ' '.join('%.2f' % x for x in R['rel'])))
def wcorr(a, b): return (w * a * b).sum() / np.sqrt((w * a * a).sum() * (w * b * b).sum())
def jac(A, ka, B, kb): sa = A['lab'] == ka; sb = B['lab'] == kb; return (sa & sb).sum() / float((sa | sb).sum())
spec = []
for j in range(NAX):
    C = {nm: (R['d'] * np.sqrt(w)[None]) @ V[j] for nm, R in RUNS.items()}; sd = np.sqrt(sum((R['pn'] / nr * C[nm] ** 2).sum() for nm, R in RUNS.items()))
    pole = {nm: (int(np.argmin(C[nm])), int(np.argmax(C[nm]))) for nm in RUNS}
    print('== axis %d (%.1f %% of the pooled leaf variance): each run\'s leaves in pooled sd units; poles = lowest / highest' % (j + 1, 100 * sh[j]))
    for nm, R in RUNS.items():
        o = np.argsort(C[nm]); lo, hi = pole[nm]
        print('  %-7s %s | poles s%02d (%.2f, %.3f of particles) / s%02d (%.2f, %.3f)' % (nm, ' '.join('s%02d:%+.2f' % (R['ks'][k], C[nm][k] / sd) for k in o), R['ks'][lo], C[nm][lo] / sd,
              R['pn'][lo], R['ks'][hi], C[nm][hi] / sd, R['pn'][hi]))
        for tag, k in (('lo', lo), ('hi', hi)): spec.append('state_%s_ax%d_%s %s %s' % (nm, j + 1, tag, R['files'][k][0], R['files'][k][1]))
    for tag, t in (('low pole', 0), ('high pole', 1)):
        for a, b in itertools.combinations(RUNS, 2):
            A, B = RUNS[a], RUNS[b]; ka, kb = pole[a][t], pole[b][t]
            rho = wcorr(A['d'][ka], B['d'][kb]) / np.sqrt(max(A['rel'][ka] * B['rel'][kb], 1e-6)); dl = (w * (A['F'][ka] - B['F'][kb]) ** 2).sum() / (A['nu'][ka] + B['nu'][kb])
            alt = [k for k in np.argsort(np.abs(C[b] - C[a][ka])) if k != kb][0]   # run B's closest non-pole leaf on this axis (control)
            rho_c = wcorr(A['d'][ka], B['d'][alt]) / np.sqrt(max(A['rel'][ka] * B['rel'][alt], 1e-6)); dl_c = (w * (A['F'][ka] - B['F'][alt]) ** 2).sum() / (A['nu'][ka] + B['nu'][alt])
            print('  %-9s %s s%02d vs %s s%02d: rho %.3f delta %.2f Jaccard %.2f | control (B s%02d) rho %.3f delta %.2f' % (tag, a, A['ks'][ka], b, B['ks'][kb], rho, dl, jac(A, A['ks'][ka], B, B['ks'][kb]),
                  B['ks'][alt], rho_c, dl_c))
    for nm, R in RUNS.items():
        lo, hi = pole[nm]; print('  within %-7s low vs high pole: rho %.3f delta %.1f' % (nm, wcorr(R['d'][lo], R['d'][hi]) / np.sqrt(max(R['rel'][lo] * R['rel'][hi], 1e-6)),
                                                                                    (w * (R['F'][lo] - R['F'][hi]) ** 2).sum() / (R['nu'][lo] + R['nu'][hi])))
if SPEC: open(SPEC, 'w').write('\n'.join(spec) + '\n'); print('spec of the extreme leaves: %s' % SPEC)
if FIG and NAX >= 2:   # every run's leaves in the pooled stable plane (sd units), marker area = leaf fraction, poles outlined
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    P2 = {nm: np.array([(R['d'] * np.sqrt(w)[None]) @ V[j] for j in range(2)]).T for nm, R in RUNS.items()}
    sds = [np.sqrt(sum((R['pn'] / nr * P2[nm][:, j] ** 2).sum() for nm, R in RUNS.items())) for j in range(2)]
    f, ax = plt.subplots(figsize=(5.2, 4.6)); col = dict(zip(RUNS, ['#2a6fdb', '#d9822b', '#3a9d5d', '#8e5bd0']))
    for nm, R in RUNS.items():
        xy = P2[nm] / np.array(sds); ax.scatter(xy[:, 0], xy[:, 1], s=2500 * R['pn'], color=col[nm], alpha=0.55, label=nm, edgecolors='none')
        for j in range(2):
            for k in (int(np.argmin(P2[nm][:, j])), int(np.argmax(P2[nm][:, j]))): ax.scatter(xy[k, 0], xy[k, 1], s=2500 * R['pn'][k], facecolors='none', edgecolors=col[nm], lw=1.4)
    ax.axhline(0, color='0.85', lw=0.8); ax.axvline(0, color='0.85', lw=0.8); ax.set_xlabel('stable axis 1 (pooled sd)'); ax.set_ylabel('stable axis 2 (pooled sd)')
    ax.legend(frameon=False, fontsize=8, markerscale=0.4); ax.set_title('r1 leaves of four independent runs (lp %.0f, node mask); poles outlined' % LP, fontsize=8); f.tight_layout(); f.savefig(FIG, dpi=150)
    print('figure %s' % FIG)

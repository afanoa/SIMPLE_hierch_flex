#!/usr/bin/env python3
"""Site discovery A on a state tree (rflex13 / rflex13s run dir): every finished node's state sets (the cheap children along the view-
centred flex PC1 = the node falloff mask's source, and the node's refined states) give a cross-half variance map v = sum_k w_k (E_k - E)
(O_k - O) of the lp-LP, least-squares scaled half maps (arm_mask2.py); its core (> FRAC of the maximum) is split into connected
components, each a CANDIDATE site with a falloff-prior mask (core dilated DIL px, weight exp(-d^2 / 2 FALL^2) beyond: arm_mask2.py).
Per candidate: volume, share of the node's reproducible variance, and the cross-half reproducibility of the node's state deviations
inside the mask, r = sum_k w_k <dE_k, dO_k>_M / sqrt(sum_k w_k |dE_k|^2_M sum_k w_k |dO_k|^2_M). FACTOR SPLIT (default): the core voxels' deviation profiles over the node's states (sqrt(w_k) x the half-averaged deviation) are clustered
by sign-invariant spherical 2-means; two clusters whose profile directions have |cos| < FACTOR_COS (each >= FACTOR_SHARE of the variance) vary
independently across the states and become two sites; one direction (incl. its sign-flipped part: an arm leaving one position for another)
stays one site. SIGNED SPLIT (--signsplit > 0, off by default): when the seed state's deviation inside a component has both signs (each >= SIGNSPLIT of its energy: density gained in one
place, lost in another, e.g. an arm moving from an open to a closed position) the component is replaced by its two parts, each an occupancy
site with its own template. The SEED of a candidate = the state of its set with the largest cross-half deviation energy inside the mask (>= SEEDMIN particles).
Candidates are de-duplicated across
nodes and sources by mask overlap (Dice of the > 0.5 regions >= DICE: one site; the representative = the most reproducible member) and
ranked by r. Writes OUT/site_NN.mrc + OUT/sites.txt (site, representative node:source:component, members, r, volume) and, with --truth
masks, recall (true site covered >= 0.8 by a discovered mask), precision and the radial coverage of every true site.
usage: discover_sites.py --tree DIR --out DIR [--truth A.mrc,B.mrc,C.mrc --frac 0.2 --dil 3 --fall 6 --lp 15 --minvox 200 --dice 0.3]"""
import os, glob, argparse, numpy as np
from scipy import ndimage
ap = argparse.ArgumentParser(); ap.add_argument('--tree', required=True); ap.add_argument('--out', required=True); ap.add_argument('--truth', default='')
ap.add_argument('--frac', type=float, default=0.2); ap.add_argument('--dil', type=int, default=3); ap.add_argument('--fall', type=float, default=6.); ap.add_argument('--lp', type=float, default=15.)
ap.add_argument('--sig', type=float, default=1.5); ap.add_argument('--minvox', type=int, default=200); ap.add_argument('--dice', type=float, default=0.3); ap.add_argument('--mskdiam', type=float, default=220.)
ap.add_argument('--box', type=int, default=160); ap.add_argument('--smpd', type=float, default=1.464); ap.add_argument('--minshare', type=float, default=0.05); ap.add_argument('--seedmin', type=int, default=1000); ap.add_argument('--all', type=int, default=0, help='1: also unfinished nodes (whatever state sets exist)')
ap.add_argument('--signsplit', type=float, default=0., help='split a component into its gained / lost parts when the smaller holds >= this share of the seed deviation energy (0 = off)')
ap.add_argument('--extmasks', default='', help='comma list of given site masks (e.g. the truth): one site each, seeded from the tree\'s node states (classification without discovery)'); ap.add_argument('--peel', type=int, default=0, help='rounds of peeling: blank found components, re-threshold the residual variance at FRAC of its own peak'); ap.add_argument('--peel_snr', type=float, default=3., help='stop peeling when the residual peak is below this multiple of the null (99.9 %% of the negative tail)'); ap.add_argument('--maxmol', type=float, default=1., help='drop candidates whose mask (> 0.5) covers more than this share of the molecule'); ap.add_argument('--allcands', type=int, default=0, help='1: write every candidate mask (cand_*.mrc) and candidates.txt'); ap.add_argument('--nodemasks', type=int, default=1, help='1: every finished node\'s own falloff mask (rsmask.mrc) is a candidate too'); ap.add_argument('--factor_share', type=float, default=0.1, help='each factor must hold this share of the component variance'); ap.add_argument('--factor_cos', type=float, default=0.8, help='split a component into two voxel groups whose deviation profiles over the states have |cos| below this (1 = off)'); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
def wr(p, x, smpd):
    x = np.ascontiguousarray(x, np.float32); h = np.zeros(256, np.int32); h[0:3] = x.shape; h[3] = 2; h[7:10] = x.shape; hf = h.view(np.float32); hf[10:13] = np.array(x.shape) * smpd; hf[13:16] = 90.
    h[16:19] = [1, 2, 3]; hf[19], hf[20], hf[21] = x.min(), x.max(), x.mean(); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = 0x4444
    with open(p, 'wb') as f: f.write(h.tobytes()); f.write(x.tobytes())
n = a.box; s = a.smpd; f = np.fft.fftfreq(n, d=s); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n, d=s), indexing='ij'); k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)
lpf = np.exp(-(k * a.lp) ** 2 * np.log(2)); lpr = lambda x: np.fft.irfftn(np.fft.rfftn(x) * lpf, s=x.shape)
g = np.arange(n) - n // 2; R = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2); sph = R <= a.mskdiam / 2 / s
gw = np.minimum(np.arange(n), n - np.arange(n)); Rw = np.sqrt(gw[:, None, None] ** 2 + gw[None, :, None] ** 2 + gw[None, None, :] ** 2)
gk = np.fft.rfftn(np.exp(-Rw ** 2 / (2 * a.sig ** 2)) / np.exp(-Rw ** 2 / (2 * a.sig ** 2)).sum()); smooth = lambda x: np.fft.irfftn(np.fft.rfftn(x) * gk, s=x.shape)
def pad(v):   # box-m map -> box n (Fourier pad, same physical size)
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); P = np.zeros((n, n, n), complex); lo = n // 2 - m // 2; P[lo:lo + m, lo:lo + m, lo:lo + m] = F; return np.real(np.fft.ifftn(np.fft.ifftshift(P))) * (n / m) ** 3
def state_set(listing):   # [(even, odd, weight)] -> lp, scaled half maps
    w = np.array([t[2] for t in listing], float); w /= w.sum(); E = np.array([lpr(pad(rd(t[0]))) for t in listing]); O = np.array([lpr(pad(rd(t[1]))) for t in listing])
    Em, Om = np.tensordot(w, E, 1), np.tensordot(w, O, 1); sup = ((Em + Om) / 2 > 0.1 * ((Em + Om) / 2)[sph].max()) & sph
    for j in range(len(w)): E[j] *= (E[j][sup] @ Em[sup]) / max(E[j][sup] @ E[j][sup], 1e-30); O[j] *= (O[j][sup] @ Om[sup]) / max(O[j][sup] @ O[j][sup], 1e-30)
    Em, Om = np.tensordot(w, E, 1), np.tensordot(w, O, 1); return w, E - Em, O - Om, sup
def repro(w, dE, dO, M):
    num = sum(wk * ((e * M) * (o * M)).sum() for wk, e, o in zip(w, dE, dO)); de = sum(wk * ((e * M) ** 2).sum() for wk, e in zip(w, dE)); do = sum(wk * ((o * M) ** 2).sum() for wk, o in zip(w, dO))
    return num / max(np.sqrt(de * do), 1e-30)
cands = []; EXT = [np.clip(rd(p), 0, 1) * sph for p in a.extmasks.split(',')] if a.extmasks else []
for D in sorted(glob.glob(a.tree + '/node_*')):
    u = D.split('/node_')[-1]
    if not os.path.exists(D + '/node.done') and not a.all: continue
    sets = []
    CR = sorted(glob.glob(D + '/cheap/rec/*_reconstruct3D'))
    if CR and os.path.exists(D + '/cheap/children_1.txt'):
        ch = np.loadtxt(D + '/cheap/children_1.txt', dtype=int); L = [(CR[-1] + '/recvol_state%02d_even_unfil.mrc' % c, CR[-1] + '/recvol_state%02d_odd_unfil.mrc' % c, (ch == c).sum(), c) for c in range(1, 11)]
        L = [t for t in L if os.path.exists(t[0]) and t[2] > 0]
        if len(L) >= 4: sets.append(('cheap', L, D + '/cheap/children_1.txt'))
    RD = sorted(glob.glob(D + '/fstageB/*_refine3D'), key=os.path.getmtime)
    if RD and os.path.exists(D + '/refined_f.txt'):
        st = np.loadtxt(D + '/refined_f.txt', dtype=int); L = [(RD[-1] + '/recvol_state%02d_even_unfil.mrc' % c, RD[-1] + '/recvol_state%02d_odd_unfil.mrc' % c, (st == c).sum(), int(c)) for c in np.unique(st[st > 0])]
        L = [t for t in L if os.path.exists(t[0]) and t[2] > 0]
        if len(L) >= 2: sets.append(('states', L, D + '/refined_f.txt'))
    for src, L, labf in sets:
        w, dE, dO, sup = state_set(L); v = smooth(sum(wk * e * o for wk, e, o in zip(w, dE, dO))) * sph; core = v > a.frac * v.max()
        okk = np.array([t[2] >= a.seedmin for t in L])
        for xi, M in enumerate(EXT):   # external masks: every state set offers a seed (the state expressing the mask most)
            ek = np.array([((e * M) * (o * M)).sum() for e, o in zip(dE, dO)]); ks = int(np.argmax(np.where(okk, ek, -np.inf))); bm = M > 0.5
            cands.append(dict(node=u, src='ext%d%s' % (xi + 1, src), comp='0', M=M, vox=int(bm.sum()), mol=float((bm & sup).sum() / sup.sum()), share=1., r=float(repro(w, dE, dO, M)), cen=np.zeros(3), rad=0.,
                              nstates=len(L), nmemb=int(sum(t[2] for t in L)), labf=labf, seed=L[ks][3], seedn=int(L[ks][2]), seedshare=float(ek[ks] / max(ek[ek > 0].sum(), 1e-30)), ek=float(ek[ks])))
            print('external mask %d at %s:%s: seed state %d (%d ptcls), deviation energy %.3g (%.0f %% of the set), repro %.3f' % (xi + 1, u, src, L[ks][3], L[ks][2], ek[ks], 100 * cands[-1]['seedshare'], cands[-1]['r']), flush=True)
        if EXT: continue
        if a.nodemasks and src == 'states' and os.path.exists(D + '/rsmask.mrc'):   # the node's own falloff mask as a candidate (discovery A in its plain form)
            M = np.clip(rd(D + '/rsmask.mrc'), 0, 1) * sph; ek = np.array([((e * M) * (o * M)).sum() for e, o in zip(dE, dO)]); ks = int(np.argmax(np.where(okk, ek, -np.inf))); bm = M > 0.5
            cands.append(dict(node=u, src='nodemask', comp='0', M=M, vox=int(bm.sum()), mol=float((bm & sup).sum() / sup.sum()), share=1., r=float(repro(w, dE, dO, M)), cen=np.array([ix.mean() for ix in np.nonzero(bm)]) if bm.any() else np.zeros(3),
                              rad=0., nstates=len(L), nmemb=int(sum(t[2] for t in L)), labf=labf, seed=L[ks][3], seedn=int(L[ks][2]), seedshare=float(ek[ks] / max(ek[ek > 0].sum(), 1e-30))))
            cands[-1]['rad'] = float(np.linalg.norm(cands[-1]['cen'] - n / 2) * s)
            print('candidate %s:nodemask:0  %d states | mask > 0.5 %.1f %% of the molecule, repro %.3f, seed state %d' % (u, len(L), 100 * cands[-1]['mol'], cands[-1]['r'], cands[-1]['seed']), flush=True)
        lab, nl = ndimage.label(core); tot = v[core].sum(); ncore = nl
        if a.peel > 0:   # PEEL: blank the found components (dilated) and re-threshold the residual variance against its own peak, while that
            vr = v.copy(); floor = np.quantile(-v[sph], 0.999); taken = core.copy()   # peak stays above the null (negative tail) by PEEL_SNR
            for _ in range(a.peel):
                vr[ndimage.binary_dilation(taken, iterations=a.dil + 2)] = 0.
                if vr.max() < a.peel_snr * max(floor, 1e-30): break
                c2 = (vr > a.frac * vr.max()) & ~taken; l2, n2 = ndimage.label(c2)
                for q in range(1, n2 + 1):
                    cq = l2 == q
                    if cq.sum() >= a.minvox: nl += 1; lab[cq] = nl
                print('   peel at %s:%s: residual peak %.3g (%.1f x the null), %d new component(s)' % (u, src, vr.max(), vr.max() / max(floor, 1e-30), n2))
                taken |= c2; tot = max(tot, 1e-30)
        def add(cm, tag, prof=None):
            M = np.exp(-ndimage.distance_transform_edt(~ndimage.binary_dilation(cm, iterations=a.dil)) ** 2 / (2 * a.fall ** 2)) * sph
            cen = np.array([ix.mean() for ix in np.nonzero(cm)]); r_ = repro(w, dE, dO, M)
            ek = np.array([((e * M) * (o * M)).sum() for e, o in zip(dE, dO)]); ks = int(np.argmax(np.where(okk, ek, -np.inf)))   # seed: the state expressing the site most
            if prof is not None: ks = int(np.argmax(np.where(okk, np.abs(prof), -np.inf)))   # a factor: the state loading most on its own profile
            cands.append(dict(node=u, src=src, comp=tag, M=M, vox=int((M > 0.5).sum()), mol=float(((M > 0.5) & sup).sum() / sup.sum()), share=float(v[cm].sum() / tot), r=float(r_),
                              cen=cen, rad=float(np.linalg.norm(cen - n / 2) * s), nstates=len(L), nmemb=int(sum(t[2] for t in L)), labf=labf, seed=L[ks][3], seedn=int(L[ks][2]), seedshare=float(ek[ks] / max(ek[ek > 0].sum(), 1e-30))))
            print('candidate %s:%s:%s  %d states %d ptcls | core %d voxels, %.0f %% of the node variance, mask > 0.5 %.1f %% of the molecule, centre %.0f A from the box centre, repro %.3f, seed state %d' %
                  (u, src, tag, len(L), cands[-1]['nmemb'], cm.sum(), 100 * cands[-1]['share'], 100 * cands[-1]['mol'], cands[-1]['rad'], r_, cands[-1]['seed']), flush=True); return ks
        for cidx in range(1, nl + 1):
            cm = lab == cidx
            if cm.sum() < a.minvox or (v[cm].sum() / tot < a.minshare and cidx <= ncore): continue
            if a.factor_cos < 1 and cm.sum() >= 2 * a.minvox:   # FACTOR SPLIT: voxel profiles over the states, sign-invariant 2-means
                idx = np.nonzero(cm); Pf = (np.sqrt(w)[:, None] * 0.5 * (dE[:, idx[0], idx[1], idx[2]] + dO[:, idx[0], idx[1], idx[2]])).T; Pn = Pf / np.maximum(np.linalg.norm(Pf, axis=1, keepdims=True), 1e-30)
                U_, S_, Vt = np.linalg.svd(Pn, full_matrices=False); Cc = Vt[:2].copy()
                for _ in range(20):
                    asg = np.argmax(np.abs(Pn @ Cc.T), 1)
                    for j in range(2):
                        if (asg == j).sum() > 2: Cc[j] = np.linalg.svd(Pn[asg == j], full_matrices=False)[2][0]
                cos = abs(float(Cc[0] @ Cc[1])); vv = v[idx]; sh_ = [vv[asg == j].sum() / vv.sum() for j in range(2)]
                if cos < a.factor_cos and min(sh_) >= a.factor_share and min((asg == j).sum() for j in range(2)) >= a.minvox:
                    print('   component %s:%s:%d splits into two factors: |cos| of their state profiles %.2f, variance shares %.2f / %.2f' % (u, src, cidx, cos, sh_[0], sh_[1]))
                    for j in range(2):
                        b = np.zeros_like(cm); b[tuple(ix[asg == j] for ix in idx)] = True; lb, nb = ndimage.label(b); sz = np.bincount(lb.ravel())[1:]
                        add(np.isin(lb, 1 + np.where(sz >= max(20, a.minvox // 4))[0]), '%df%d' % (cidx, j + 1), prof=Cc[j] / np.sqrt(w))
                    continue
                print('   component %s:%s:%d is one factor (|cos| %.2f, shares %.2f / %.2f)' % (u, src, cidx, cos, sh_[0], sh_[1]))
            ek = np.array([((e * cm) * (o * cm)).sum() for e, o in zip(dE, dO)]); ks = int(np.argmax(np.where(okk, ek, -np.inf))); dk = 0.5 * (dE[ks] + dO[ks])
            ep, en = (dk[cm & (dk > 0)] ** 2).sum(), (dk[cm & (dk < 0)] ** 2).sum(); parts = []
            if a.signsplit > 0 and min(ep, en) / max(ep + en, 1e-30) >= a.signsplit:   # density gained AND lost: two occupancy sites
                for sg, b in (('p', cm & (dk > 0)), ('n', cm & (dk < 0))):
                    lb, nb = ndimage.label(b); sz = np.bincount(lb.ravel())[1:]
                    if nb and sz.max() >= a.minvox: parts.append((sg, lb == (np.argmax(sz) + 1)))
            if len(parts) == 2:
                print('   component %s:%s:%d split by the sign of the seed state %d deviation (energy + %.0f %% / - %.0f %%)' % (u, src, cidx, L[ks][3], 100 * ep / (ep + en), 100 * en / (ep + en)))
                for sg, b in parts: add(b, '%d%s' % (cidx, sg))
            else: add(cm, '%d' % cidx)
if a.allcands:   # every candidate's mask + table (same columns as sites.txt) for hand-picked site sets
    with open(a.out + '/candidates.txt', 'w') as fo:
        for i, c0 in enumerate(cands, 1):
            fn = '%s/cand_%s_%s_%s.mrc' % (a.out, c0['node'], c0['src'], c0['comp']); wr(fn, c0['M'], s)
            fo.write('%d %s:%s:%s %.3f %d %.3f %.1f %s %d %d %s\n' % (i, c0['node'], c0['src'], c0['comp'], c0['r'], c0['vox'], c0['mol'], c0['rad'], c0['labf'], c0['seed'], c0['seedn'], fn))
if EXT:   # one site per external mask: the state set whose best state carries the most deviation energy inside it
    keep = [max((c for c in cands if c['src'].startswith('ext%d' % (xi + 1))), key=lambda c: c['ek']) for xi in range(len(EXT))]
    cands = keep
    for c in cands: print('external mask site %s: seed state %d of %s (%d ptcls)' % (c['src'], c['seed'], c['labf'], c['seedn']))
if a.maxmol < 1:   # a site is local: drop candidates covering more than MAXMOL of the molecule (global / junk-driven variation)
    nb = len(cands); cands = [c for c in cands if c['mol'] <= a.maxmol]; print('dropped %d of %d candidates covering > %.0f %% of the molecule' % (nb - len(cands), nb, 100 * a.maxmol))
# ---- de-duplicate (greedy by reproducibility)
order = sorted(range(len(cands)), key=lambda i: (cands[i]['src'] == 'nodemask', -cands[i]['r'])); sites = []   # a node's whole mask is the fallback representative
def dice(m1, m2): b1, b2 = m1 > 0.5, m2 > 0.5; return 2. * (b1 & b2).sum() / max(b1.sum() + b2.sum(), 1)
for i in order:
    hit = [j for j, S in enumerate(sites) if dice(cands[i]['M'], cands[S[0]]['M']) >= a.dice]
    if hit: sites[hit[0]].append(i)
    else: sites.append([i])
with open(a.out + '/sites.txt', 'w') as fo:
    fo.write('# site representative(node:source:component) repro volume_vox molecule_share centre_radius_A seed_label_file seed_id seed_particles members\n')
    for j, S in enumerate(sites, 1):
        c0 = cands[S[0]]; wr('%s/site_%02d.mrc' % (a.out, j), c0['M'], s)
        fo.write('%d %s:%s:%s %.3f %d %.3f %.1f %s %d %d %s\n' % (j, c0['node'], c0['src'], c0['comp'], c0['r'], c0['vox'], c0['mol'], c0['rad'], c0['labf'], c0['seed'], c0['seedn'], ' '.join('%s:%s:%s' % (cands[i]['node'], cands[i]['src'], cands[i]['comp']) for i in S)))
        print('SITE %d: %s:%s:%s repro %.3f, mask > 0.5 %d voxels (%.1f %% of the molecule), centre %.0f A; seed = state %d (%d ptcls, %.0f %% of the site deviation energy); %d member candidates (%s)' % (j, c0['node'], c0['src'], c0['comp'], c0['r'], c0['vox'], 100 * c0['mol'], c0['rad'], c0['seed'], c0['seedn'], 100 * c0['seedshare'], len(S), ' '.join('%s:%s' % (cands[i]['node'], cands[i]['src']) for i in S[:8])))
# ---- truth
if a.truth:
    TM = [rd(p) for p in a.truth.split(',')]; tn = [os.path.basename(p).replace('.mrc', '') for p in a.truth.split(',')]
    SM = [rd('%s/site_%02d.mrc' % (a.out, j)) for j in range(1, len(sites) + 1)]
    for t, nm in zip(TM, tn):
        tb = t > 0.5; cov = [((m > 0.5) & tb).sum() / tb.sum() for m in SM]; j = int(np.argmax(cov)) if cov else -1
        cc_ = sorted([(((cands[i]['M'] > 0.5) & tb).sum() / tb.sum(), ((cands[i]['M'] > 0.5) & (t > 0.05)).sum() / max((cands[i]['M'] > 0.5).sum(), 1), i) for i in range(len(cands))], reverse=True); cands_cov = cc_[0][0] if cc_ else 0.
        print('   %s: best candidates (coverage, share of the candidate inside): %s' % (nm, '; '.join('%s:%s:%s %.2f / %.2f' % (cands[i]['node'], cands[i]['src'], cands[i]['comp'], c1, c2) for c1, c2, i in cc_[:5])))
        dd_ = sorted([(dice(cands[i]['M'], t), i) for i in range(len(cands))], reverse=True)
        print('   %s: best candidates by Dice with the truth (> 0.5): %s' % (nm, '; '.join('%s:%s:%s %.2f' % (cands[i]['node'], cands[i]['src'], cands[i]['comp'], d_) for d_, i in dd_[:4])))
        print('TRUTH %s (%d voxels > 0.5): best site %d covers %.2f (recall %s), its mask is %.0f %% inside the truth (> 0.05); best of all candidates %.2f' %
              (nm, tb.sum(), j + 1, cov[j] if cov else 0, 'yes' if cov and cov[j] >= 0.8 else 'no', 100 * ((SM[j] > 0.5) & (t > 0.05)).sum() / max((SM[j] > 0.5).sum(), 1) if cov else 0, cands_cov))
    anyt = sum(TM) > 0.05; prec = [((m > 0.5) & anyt).sum() / max((m > 0.5).sum(), 1) for m in SM]
    print('PRECISION: %d of %d discovered sites lie >= 50 %% inside a true site (fraction inside per site: %s)' % (sum(p >= 0.5 for p in prec), len(SM), ' '.join('%.2f' % p for p in prec)))
    edges = np.arange(0, 125, 10.); rb = np.digitize(R * s, edges) - 1
    for t, nm in zip(TM, tn):
        j = int(np.argmax([((m > 0.5) & (t > 0.5)).sum() for m in SM])); line = []
        for i in range(len(edges) - 1):
            m_ = (t > 0.5) & (rb == i)
            if m_.sum() >= 20: line.append('%.0f-%.0f A %.2f' % (edges[i], edges[i + 1], SM[j][m_].mean()))
        print('RADIAL %s by site %d (mean mask value on the true site per radius bin): %s' % (nm, j + 1, ', '.join(line)))

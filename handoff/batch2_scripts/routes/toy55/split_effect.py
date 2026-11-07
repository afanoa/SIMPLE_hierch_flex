#!/usr/bin/env python3
"""usage: split_effect.py TREE_DIR DUMP_NPZ OUT [--fa MAPS_TABLE] [--dbref QS3_FILE] [--res 15]
Effect size of every split of a hetero_tree run. The split groups come from the run's tree.log ('node X: split group G (N particles;
stage-B states s:n ..) -> node Y'), each group's map is node X's stage-B map of the state that dominates it (stageB/recvol_state0S
with its even / odd halves), its rows are node X's split/split_children.txt == G. Every group is compared with the largest one:
- z: the split-half correlation of D = M_g - M_ref over shells 1..RES in the soft node mask (Re<E_g - E_ref, O_g - O_ref> over the
  norms) x sqrt(N_eff); snr: |D|^2 over its half-map noise;
- share: the share of the mask where the low-passed D exceeds 3 sd of its half-map noise (how much of the molecule changes);
- loc: D less its per-shell projection on the two groups' mean map, its reproducible power over the mean map's (phaseF_rel.py's
  rel_loc); cos: D's signal cosine to the mean map (1 = a rescaled copy: brightness, dilution);
- top: the centroid of the top 1 % abs(D), A from the box centre, and its radius from the box centre; +/-: D's sign balance;
- qAUC: corr of group g against the reference group at matched defocus (dfx deciles; root stage-B corr from DUMP_NPZ), 0.5 = equal;
- with --fa (fa_maps10.py's maps_table.txt holding rows named split_X_s<S> for the stage-B maps, split_maps.py's list): ds1 (A),
  dclosed and dcore (core CC mean), group g minus the reference;
- dB: quality_screen3.py's blur of the lower-corr group in the stable region (mask voxels with abs(D) under 3 sd), from the
  UNFILTERED halves' cross-power (stage-B full maps are filtered per state): -4 x the slope of 0.5 log(S_L / S_H) against s^2 over
  15 A -> max(6 A, 2.5 px), A^2; with --dbref (quality_ref's QS3 line: the root's top vs bottom corr quintile at the same box) also
  dB_ref, ratio = dB / dB_ref and blur = 1 when ratio >= 1/2 and dB >= 3 SE (validated 2026-10-25: quality 0.51..0.64, real <= 0.18).
One SPLIT line per compared pair, then a header-commented table."""
import os, re, sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('tree'); ap.add_argument('dump'); ap.add_argument('out'); ap.add_argument('--fa'); ap.add_argument('--dbref'); ap.add_argument('--res', type=float, default=15.)
a = ap.parse_args(); T = a.tree.rstrip('/'); TL = os.path.join(os.path.dirname(T), 'tree.log')
def read_mrc(f):
    raw = open(f, 'rb').read(); h = np.frombuffer(raw[:1024], np.int32); hf = h.view(np.float32); n = int(h[0])
    return np.frombuffer(raw[1024 + int(h[23]):], np.float32, count=n ** 3).astype(np.float64).reshape(n, n, n), float(hf[10]) / n
def fresample(v, n):
    m = v.shape[0]
    if m == n: return v
    F = np.fft.fftshift(np.fft.fftn(v)); c = m // 2
    if n < m: s = c - n // 2; G = F[s:s + n, s:s + n, s:s + n]
    else: G = np.zeros((n, n, n), complex); s = n // 2 - c; G[s:s + m, s:s + m, s:s + m] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / m) ** 3
def kr(n): f = np.fft.fftfreq(n) * n; return np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
def cdot(x, y): return float((x.conj() * y).real.sum())
D = np.load(a.dump); corr, dfx = D['corr'], D['dfx']; ok = np.isfinite(dfx)
dd = np.full(len(dfx), -1); dd[ok] = np.searchsorted(np.quantile(dfx[ok], np.linspace(0.1, 0.9, 9)), dfx[ok])
def auc(x, y):
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) < 20 or len(y) < 20: return np.nan
    r = np.argsort(np.argsort(np.r_[x, y])); return (r[:len(x)].sum() - len(x) * (len(x) - 1) / 2.) / (len(x) * len(y))
FA = {}
if a.fa and os.path.exists(a.fa):
    hd = None
    for l in open(a.fa):
        w = l.split()
        if l.startswith('#'): hd = w[1:]; continue
        if hd and w: FA[w[0]] = dict(zip(hd, w))
RB = re.search(r'dBx ([+-]?[0-9.]+) A\^2 \(SE ([0-9.]+)\)', open(a.dbref).read()) if a.dbref and os.path.exists(a.dbref) else None
RB = (float(RB.group(1)), float(RB.group(2))) if RB else (np.nan, np.nan)
SP = {}
for l in open(TL, errors='replace'):
    m = re.search(r'HETERO_TREE node (\S+): split group (\d+) \((\d+) particles; stage-B states ([^)]*)\) -> (.*)$', l)
    if not m: continue
    st = {int(p.split(':')[0]): int(p.split(':')[1]) for p in m.group(4).split()}; SP.setdefault(m.group(1), {})[int(m.group(2))] = (int(m.group(3)), st, m.group(5).strip().replace(' ', '_'))
rows = []
for X in sorted(SP, key=lambda s: (len(s), s)):
    N = '%s/node_%s' % (T, X); G = SP[X]
    if not os.path.exists(N + '/split/split_children.txt'): continue
    ch = np.loadtxt(N + '/split/split_children.txt', dtype=int); ref = max(G, key=lambda g: G[g][0]); dom = {g: max(G[g][1], key=G[g][1].get) for g in G}
    if len(ch) != len(corr): c2 = np.zeros(len(corr), int); c2[np.loadtxt(N + '/rows.txt', dtype=int) - 1] = ch; ch = c2  # newer builds: per node row
    mp = {}
    for g in G:
        s = dom[g]; f = ['%s/stageB/recvol_state%02d%s.mrc' % (N, s, x) for x in ('', '_even', '_odd')]
        if all(os.path.exists(p) for p in f): mp[g] = [read_mrc(p)[0] for p in f]
    if ref not in mp: continue
    smpd = read_mrc('%s/stageB/recvol_state%02d.mrc' % (N, dom[ref]))[1]; n = mp[ref][0].shape[0]
    mk = (fresample(read_mrc(N + '/mask.mrc')[0], n) > 0.5).astype(float); sh = np.rint(kr(n)).astype(int).ravel(); smax = int(n * smpd / a.res); sel = (sh >= 1) & (sh <= smax); shs = sh[sel]
    sm = np.real(np.fft.ifftn(np.fft.fftn(mk) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2))); neff = mk.sum() * smpd ** 3 / (a.res / 2) ** 3
    lpw = lambda: np.clip((1.1 * smpd / a.res - kr(n) / n) / (0.2 * smpd / a.res), 0., 1.)
    w = 0.5 - 0.5 * np.cos(np.pi * lpw()); g_ = np.arange(n) - n // 2; ZZ, YY, XX = np.meshgrid(g_, g_, g_, indexing='ij')
    Q = {g: [np.fft.fftn(x * sm).ravel()[sel] for x in mp[g]] for g in mp}
    for g in sorted(mp):
        if g == ref: continue
        d = [Q[g][j] - Q[ref][j] for j in range(3)]; pm = [(Q[g][j] + Q[ref][j]) / 2 for j in range(3)]; rep = cdot(d[1], d[2])
        z = rep / np.sqrt(cdot(d[1], d[1]) * cdot(d[2], d[2])) * np.sqrt(neff); nd = ((np.abs(Q[g][1] - Q[g][2]) ** 2 + np.abs(Q[ref][1] - Q[ref][2]) ** 2) / 4).sum()
        snr = (np.abs(d[0]) ** 2).sum() / nd; sigm = cdot(pm[1], pm[2]); cs = cdot(d[1], pm[2]) / np.sqrt(max(rep, 1e-30) * max(sigm, 1e-30)) if rep > 0 else np.nan
        num = np.bincount(shs, (d[0].conj() * pm[0]).real, smax + 1); den = np.bincount(shs, np.abs(pm[0]) ** 2, smax + 1); c = (num / np.maximum(den, 1e-30))[shs]
        loc = cdot(d[1] - c * pm[1], d[2] - c * pm[2]) / sigm
        Dr = np.real(np.fft.ifftn((np.fft.fftn(mp[g][0]) - np.fft.fftn(mp[ref][0])) * w)); Nr = np.real(np.fft.ifftn(((np.fft.fftn(mp[g][1]) - np.fft.fftn(mp[g][2])) - (np.fft.fftn(mp[ref][1]) - np.fft.fftn(mp[ref][2]))) * w)) / 2
        m_ = mk > 0.5; share = (np.abs(Dr[m_]) > 3 * Nr[m_].std()).mean(); top = m_ & (np.abs(Dr) >= np.percentile(np.abs(Dr[m_]), 99))
        cx, cy, cz = XX[top].mean() * smpd, YY[top].mean() * smpd, ZZ[top].mean() * smpd; pmr = (Dr[m_] ** 2 * (Dr[m_] > 0)).sum() / max((Dr[m_] ** 2 * (Dr[m_] < 0)).sum(), 1e-30)
        q = np.nanmean([auc(corr[(ch == g) & (dd == j)], corr[(ch == ref) & (dd == j)]) for j in range(10)])
        Lg = np.nanmedian(corr[ch == g]) < np.nanmedian(corr[ch == ref]); stab = m_ & (np.abs(Dr) <= 3 * Nr[m_].std())
        if stab.sum() < 0.1 * m_.sum(): stab = m_ & (np.abs(Dr) <= np.percentile(np.abs(Dr[m_]), 30))
        ss = np.real(np.fft.ifftn(np.fft.fftn(stab.astype(float)) * np.exp(-0.5 * (2. * 2 * np.pi * kr(n) / n) ** 2))); k3 = np.rint(kr(n)).astype(int)
        uf = {h: ['%s/stageB/recvol_state%02d_%s_unfil.mrc' % (N, dom[h], x) for x in ('even', 'odd')] for h in (g, ref)}; dB = dBse = np.nan
        if all(os.path.exists(p) for h in uf for p in uf[h]):
            gl, gh = (g, ref) if Lg else (ref, g); EL, OL, EH, OH = (np.fft.fftn(read_mrc(p)[0] * ss) for p in uf[gl] + uf[gh]); xs, ys = [], []
            for k in range(int(np.ceil(n * smpd / 15.)), int(n * smpd / max(6., 2.5 * smpd)) + 1):
                q_ = k3 == k; SL, SH = (EL[q_] * OL[q_].conj()).real.sum(), (EH[q_] * OH[q_].conj()).real.sum()
                if SL > 0 and SH > 0: xs.append((k / (n * smpd)) ** 2); ys.append(0.5 * np.log(SL / SH))
            if len(xs) > 3:
                xs, ys = np.array(xs), np.array(ys); c_ = np.polyfit(xs, ys, 1); r_ = ys - np.polyval(c_, xs); dB = -4 * c_[0]
                dBse = 4 * np.sqrt(r_ @ r_ / (len(xs) - 2) / ((xs - xs.mean()) ** 2).sum())
        fa_ = [np.nan] * 3; ka, kb = 'split_%s_s%d' % (X, dom[g]), 'split_%s_s%d' % (X, dom[ref])
        if ka in FA and kb in FA: fa_ = [float(FA[ka][k]) - float(FA[kb][k]) for k in ('s1_A', 'closed_arm', 'core_cc_mean')]
        bl = int(bool(np.isfinite(dB) and RB[0] >= 3 * RB[1] and dB >= RB[0] / 2 and dB >= 3 * dBse))
        rows.append((X, g, G[g][0], G[g][2], ref, G[ref][0], G[ref][2], z, snr, share, loc, cs, pmr, cx, cy, cz, np.sqrt(cx * cx + cy * cy + cz * cz), q, *fa_, dB, g if Lg else ref, dBse, RB[0], dB / RB[0], bl))
        print('SPLIT %s group %d (%d, %s) vs %d (%d, %s) | z %.1f snr %.1f share %.3f loc %.4f cos %+.3f +/- %.2f | top (%.0f, %.0f, %.0f) r %.0f A | qAUC %.3f | ds1 %+.2f dclosed %+.3f dcore %+.3f'
              ' | dB %+.1f A^2 (low-corr group %d, SE %.1f)' % rows[-1][:24] + ' | dB / dB_ref %.2f%s' % (rows[-1][25], ' BLUR' if bl else '')); sys.stdout.flush()
with open(a.out, 'w') as f:
    f.write('# node group n child ref n_ref child_ref z snr share loc cos pm top_x top_y top_z top_r qAUC ds1 dclosed dcore dB low_corr_group dB_se dB_ref ratio blur\n')
    for r in rows: f.write(' '.join(str(x) if isinstance(x, (str, int, np.integer)) else '%.4g' % x for x in r) + '\n')

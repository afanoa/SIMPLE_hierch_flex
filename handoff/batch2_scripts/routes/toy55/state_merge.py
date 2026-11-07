#!/usr/bin/env python3
"""Bottom-up organisation of states by gold-standard distinguishability.
Every state has independent even/odd half-maps. For states a, b the cross-half difference statistic
    num(a,b) = Re < M(E_a - E_b), M(O_a - O_b) >   over the Fourier band [hp, lp] inside the envelope M
is an unbiased estimate of the true difference energy (the halves' noise is independent), zero when a and b are the same
structure; z = num / sd with sd from the half-map noise per shell (|E - O|^2 / 2 per coefficient) inflated by the envelope's volume
fraction for the mask-induced correlation. Agglomerative: merge the least distinguishable pair (smallest z), the merged node's
half-maps the particle-weighted averages (linear, so the pair statistics update on the Gram matrix of the band-limited halves),
until one node is left. The tree is the merge order; the cut at the first merge with z >= ZTHR gives the distinguishable states.
usage: state_merge.py OUT (--rdir REFINE3D_STATES_DIR --labels per-row-states.txt | --list name/even/odd/count file)
       [--smpd 1.464 --hp 40 --lp 10 --zthr 5 --unfil --gt gt_particles.txt]"""
import os, sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--rdir', default=''); ap.add_argument('--labels', default=''); ap.add_argument('--list', default='')
ap.add_argument('--smpd', type=float, default=0.); ap.add_argument('--hp', type=float, default=40.); ap.add_argument('--lp', type=float, default=10.)
ap.add_argument('--zthr', type=float, default=5.); ap.add_argument('--unfil', action='store_true'); ap.add_argument('--by', default='d2', choices=['d2', 'z']); ap.add_argument('--ncut', default='12,8,6,4'); ap.add_argument('--arm', default=''); ap.add_argument('--link', default='average', choices=['centroid', 'average', 'complete']); ap.add_argument('--budget', type=float, default=0.); ap.add_argument('--gt', default=''); ap.add_argument('--mask', default='')
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n).astype(np.float64)
lab = None; names, ev, od, cnt = [], [], [], []
if a.list:
    for l in open(a.list):
        if l.strip() and not l.startswith('#'): w = l.split(); names.append(w[0]); ev.append(w[1]); od.append(w[2]); cnt.append(float(w[3]))
else:
    lab = np.loadtxt(a.labels, dtype=int); sfx = '_unfil' if a.unfil else ''
    for s in range(1, lab.max() + 1):
        e, o = '%s/recvol_state%02d_even%s.mrc' % (a.rdir, s, sfx), '%s/recvol_state%02d_odd%s.mrc' % (a.rdir, s, sfx)
        if (lab == s).sum() > 0 and os.path.isfile(e) and os.path.isfile(o): names.append(str(s)); ev.append(e); od.append(o); cnt.append(float((lab == s).sum()))
N = len(names); cnt = np.array(cnt); E0 = rd(ev[0]); n = E0.shape[0]
if a.smpd <= 0: a.smpd = float(np.fromfile(ev[0], np.float32, 13)[10]) / n   # from the header: cell length / nx
if a.smpd <= 0: sys.exit('no pixel size in the map header: pass --smpd')
print('%d states, box %d, band %.0f - %.0f A, z threshold %.1f, %s linkage on %s' % (N, n, a.hp, a.lp, a.zthr, a.link, a.by))
# envelope: particle-weighted average of the full maps, low-passed to 15 A, thresholded at 10 % of its maximum, dilated 4 px, soft 2 px
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, np.fft.rfftfreq(n), indexing='ij'); kr = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2)   # cycles / px
if a.mask: M = np.clip(rd(a.mask), 0, 1)
else:
    avg = sum(c * (rd(e) + rd(o)) / 2 for c, e, o in zip(cnt, ev, od)) / cnt.sum()
    lpf = np.exp(-(kr * a.smpd * 15.) ** 2 * 2); sm = np.fft.irfftn(np.fft.rfftn(avg) * lpf, s=avg.shape); b = (sm > 0.1 * sm.max()).astype(float)
    r = np.sqrt(sum(np.minimum(np.arange(n), n - np.arange(n))[ix] ** 2 for ix in np.ix_(range(n), range(n), range(n))))
    ball = np.fft.rfftn((r <= 4).astype(float)); b = (np.fft.irfftn(np.fft.rfftn(b) * ball, s=b.shape) > 0.5).astype(float)
    g = np.fft.rfftn(np.exp(-r ** 2 / (2 * 2. ** 2))); g /= g.flat[0]; M = np.clip(np.fft.irfftn(np.fft.rfftn(b) * g, s=b.shape), 0, 1)
fmask = float(M.mean()); infl = M.size * float((M ** 4).sum()) / float((M ** 2).sum()) ** 2; res = np.where(kr > 0, a.smpd / np.maximum(kr, 1e-9), 1e9); band = (res <= a.hp) & (res >= a.lp)
shell = np.round(kr * n).astype(int)[band]; ush, sh_idx = np.unique(shell, return_inverse=True); nsh = np.bincount(sh_idx).astype(float)
# rfftn half-space: interior coefficients stand for two (Friedel), the kz = 0 / Nyquist planes for one
wgt = np.where((kz[band] == 0) | (np.abs(kz[band]) == 0.5), 1., 2.)
print('envelope volume fraction %.3f (correlation inflation %.1f); %d band coefficients in %d shells' % (fmask, infl, band.sum(), len(ush)))
Ev = np.zeros((N, band.sum()), complex); Ov = np.zeros_like(Ev); V = np.zeros((N, len(ush)))
for i in range(N):
    e = np.fft.rfftn(M * rd(ev[i]))[band]; o = np.fft.rfftn(M * rd(od[i]))[band]; Ev[i], Ov[i] = e, o
    V[i] = np.bincount(sh_idx, weights=wgt * np.abs(e - o) ** 2 / 2, minlength=len(ush)) / np.bincount(sh_idx, weights=wgt, minlength=len(ush))   # per-coefficient half-map noise variance per shell
W = wgt[None, :]; C = np.real((Ev * W) @ Ov.conj().T)   # C[a, b] = Re <E_a, O_b> over the band (both Friedel halves)
nw = np.bincount(sh_idx, weights=wgt, minlength=len(ush))   # coefficients per shell, both Friedel halves
def zab(a_, b_, C, V):
    num = C[a_, a_] + C[b_, b_] - C[a_, b_] - C[b_, a_]; var = (nw * (V[a_] + V[b_]) ** 2).sum() * infl   # Var Re(D1 D2*) = V1 V2 / 2 per complex pair, twice per Friedel pair; the envelope correlates ~infl neighbours
    return num / np.sqrt(var)
def d2ab(a_, b_, C):   # relative squared difference of the true maps: unbiased difference energy over the geometric mean of the unbiased signal energies
    num = C[a_, a_] + C[b_, b_] - C[a_, b_] - C[b_, a_]; return num / np.sqrt(max(C[a_, a_], 1e-30) * max(C[b_, b_], 1e-30))
# truth readout
conf = None
if a.gt and lab is not None:
    G = np.loadtxt(a.gt, comments='#'); conf = np.full(int(G[:, 0].max()) + 1, -1); conf[G[:, 0].astype(int)] = G[:, 1].astype(int); conf = conf[1:len(lab) + 1]
def cls(k): return 'junk' if k == 0 else 'grid' if k >= 400 else 'breath' if k >= 300 else 'swing' if k >= 200 else 'holo' if k >= 101 else 'apo'
arm = {}
if a.arm:   # real data: per state the closed-arm occupancy of its map (arm_states.txt: map closed open call)
    import re
    for l in open(a.arm):
        if l.strip() and not l.startswith('#'): w_ = l.split(); arm[str(int(re.findall(r'(\d+)\.mrc$', w_[0])[0]))] = (float(w_[1]), w_[3])
def coherence(groups):   # share of the leaves' variance (particle-weighted) explained by the groups: closed-arm (real) or the true swing / class (synthetic)
    out = []
    if arm:
        v = np.array([arm[m][0] if m in arm else np.nan for m in names]); n_ = cnt.copy(); ok = np.isfinite(v); mu = (n_[ok] * v[ok]).sum() / n_[ok].sum()
        tot = (n_[ok] * (v[ok] - mu) ** 2).sum(); bet = 0.
        for g in groups:
            ix = [names.index(m) for m in g if np.isfinite(v[names.index(m)])]
            if ix: bet += n_[ix].sum() * ((n_[ix] * v[ix]).sum() / n_[ix].sum() - mu) ** 2
        out.append('closed-arm eta2 %.2f' % (bet / tot))
    if conf is not None:
        gl = np.zeros(len(lab), int)
        for gi, g in enumerate(groups, 1): gl[np.isin(lab, [int(m) for m in g])] = gi
        sw = (conf >= 200) & (conf <= 230) & (gl > 0)
        if sw.sum() > 100:
            x = (conf[sw] - 200).astype(float); gg = gl[sw]; mu = x.mean(); bet = sum(((gg == k).sum() * (x[gg == k].mean() - mu) ** 2) for k in np.unique(gg))
            out.append('swing eta2 %.2f' % (bet / ((x - mu) ** 2).sum()))
        cl = np.array([cls(k) for k in conf]); on = gl > 0; pur = sum(max(((gl == k) & (cl == c)).sum() for c in set(cl[gl == k])) for k in np.unique(gl[on])) / on.sum()
        out.append('class purity %.2f' % pur)
    return ', '.join(out)
def describe(members):
    if arm:
        v = [(arm[m][0], cnt[names.index(m)], m) for m in members if m in arm]
        if not v: return ''
        mx = max(v); calls = [arm[m][1] for m in members if m in arm]
        return 'closed-arm mean %.2f max %.2f (state %s), %d closed %d intermediate %d open' % (sum(x * c for x, c, _ in v) / sum(c for _, c, _ in v), mx[0], mx[2], calls.count('closed'), calls.count('intermediate'), calls.count('open'))
    if conf is None: return ''
    rows = np.isin(lab, [int(m) for m in members]); c = conf[rows]; cl = np.array([cls(k) for k in c]); out = []
    for nm in ('apo', 'holo', 'swing', 'breath', 'grid', 'junk'):
        m = cl == nm
        if m.mean() >= 0.05:
            x = ''
            if nm == 'swing': x = ' %.1f+-%.1f A' % ((c[m] - 200).mean(), (c[m] - 200).std())
            if nm == 'breath': x = ' %.1f+-%.1f' % ((c[m] - 300).mean(), (c[m] - 300).std())
            if nm == 'grid': x = ' k %.1f j %.1f' % (((c[m] - 400) // 10).mean(), ((c[m] - 400) % 10).mean())
            out.append('%s %.0f%%%s' % (nm, 100 * m.mean(), x))
    return ', '.join(out)
# pairwise table of the leaves
Z0 = np.zeros((N, N))
for i in range(N):
    for j in range(i + 1, N): Z0[i, j] = Z0[j, i] = zab(i, j, C, V)
np.savetxt(a.out + '/leaf_z.txt', Z0, fmt='%.1f', header=' '.join(names))
D0 = np.zeros((N, N))
for i in range(N):
    for j in range(i + 1, N): D0[i, j] = D0[j, i] = d2ab(i, j, C)
np.savetxt(a.out + '/leaf_d2.txt', D0, fmt='%.5f', header=' '.join(names))
print('leaf pair z (cross-half difference energy / sd): min %.1f, median %.1f; pairs below %.1f: %d of %d' % (Z0[np.triu_indices(N, 1)].min(), np.median(Z0[np.triu_indices(N, 1)]), a.zthr, (Z0[np.triu_indices(N, 1)] < a.zthr).sum(), N * (N - 1) // 2))
# agglomerative merging on the Gram matrix
ncuts = set(int(x) for x in a.ncut.split(',') if x); alive = list(range(N)); mem = {i: [names[i]] for i in range(N)}; w = {i: cnt[i] for i in range(N)}
Dg = np.zeros((2 * N, 2 * N)); Dg[:N, :N] = D0   # cluster distances for average / complete linkage
Cg = np.zeros((2 * N, 2 * N)); Cg[:N, :N] = C; Vg = np.zeros((2 * N, len(ush))); Vg[:N] = V; steps = []; nxt = N; cut = None
phase1 = a.budget > 0   # budget phase: only a state below BUDGET particles may merge, into its nearest; then the tree over the budget states
def write_cut(tag, groups_):
    groups_ = sorted(groups_, key=lambda g: -sum(cnt[names.index(m)] for m in g)); cd = '%s/%s' % (a.out, tag); os.makedirs(cd, exist_ok=True)
    with open(cd + '/groups.txt', 'w') as fo:
        fo.write('# bottom-up state (%s linkage on %s, %s): particles, members (leaf states)\n' % (a.link, a.by, tag))
        for gi, g in enumerate(groups_, 1): fo.write('%d %d %s\n' % (gi, sum(cnt[names.index(m)] for m in g), '+'.join(g)))
    if lab is not None:
        lut = np.zeros(lab.max() + 1, int)
        for gi, g in enumerate(groups_, 1): lut[[int(m) for m in g]] = gi
        np.savetxt(cd + '/states.txt', lut[np.clip(lab, 0, None)] * (lab > 0), fmt='%d')
while len(alive) > 1:
    best = None
    for ii, p in enumerate(alive):
        for q in alive[ii + 1:]:
            if phase1 and w[p] >= a.budget and w[q] >= a.budget: continue
            z = zab(p, q, Cg, Vg); d = d2ab(p, q, Cg) if a.link == 'centroid' else Dg[p, q]; key = d if a.by == 'd2' else z
            if best is None or key < best[0]: best = (key, p, q, z, d)
    if best is None:   # every state reached the budget
        ws_ = sorted(w[k] for k in alive); gl_ = [mem[k] for k in alive]
        print('== budget %d: %d states, particles largest %d median %d smallest %d; %s' % (a.budget, len(alive), ws_[-1], ws_[len(ws_) // 2], ws_[0], coherence(gl_)))
        for k in sorted(alive, key=lambda k: -w[k]): print('   [%s] %d particles  %s' % (' '.join(mem[k]) if len(mem[k]) <= 12 else '%d leaves' % len(mem[k]), w[k], describe(mem[k])))
        write_cut('budget%d' % a.budget, gl_); phase1 = False; continue
    _, p, q, z, d = best; wp, wq = w[p], w[q]; s_ = wp + wq; r_ = nxt; nxt += 1
    for k in alive:   # the merged node's half-maps are (wp X_p + wq X_q) / s: rows and columns of the Gram matrix are linear in them
        Cg[r_, k] = (wp * Cg[p, k] + wq * Cg[q, k]) / s_; Cg[k, r_] = (wp * Cg[k, p] + wq * Cg[k, q]) / s_
    for k in alive:
        Dg[r_, k] = Dg[k, r_] = (wp * Dg[p, k] + wq * Dg[q, k]) / s_ if a.link == 'average' else max(Dg[p, k], Dg[q, k])
    Cg[r_, r_] = (wp * wp * Cg[p, p] + wp * wq * (Cg[p, q] + Cg[q, p]) + wq * wq * Cg[q, q]) / s_ ** 2
    Vg[r_] = (wp ** 2 * Vg[p] + wq ** 2 * Vg[q]) / s_ ** 2
    if cut is None and z >= a.zthr and a.by == 'z': cut = [list(mem[k]) for k in alive]
    mem[r_] = mem[p] + mem[q]; w[r_] = s_
    if len(alive) - 1 in ncuts:
        groups_ = sorted([mem[k] for k in alive if k not in (p, q)] + [mem[p] + mem[q]], key=lambda g: -sum(cnt[names.index(m)] for m in g))
        cd = '%s/cut%d' % (a.out, len(alive) - 1); os.makedirs(cd, exist_ok=True)
        with open(cd + '/groups.txt', 'w') as fo:
            fo.write('# bottom-up state (%s linkage on %s, cut %d): particles, members (leaf states)\n' % (a.link, a.by, len(alive) - 1))
            for gi, g in enumerate(groups_, 1): fo.write('%d %d %s\n' % (gi, sum(cnt[names.index(m)] for m in g), '+'.join(g)))
        if lab is not None:
            lut = np.zeros(lab.max() + 1, int)
            for gi, g in enumerate(groups_, 1): lut[[int(m) for m in g]] = gi
            np.savetxt(cd + '/states.txt', lut[np.clip(lab, 0, None)] * (lab > 0), fmt='%d')
        ws_ = sorted([w[k] for k in alive if k not in (p, q)] + [s_]); nl_ = [len(mem[k]) for k in alive if k not in (p, q)] + [len(mem[p]) + len(mem[q])]
        print('== %d states after merge %d (d2 %.4f): particles largest %d median %d smallest %d; single-leaf states %d; %s' % (len(alive) - 1, len(steps) + 1, d, ws_[-1], ws_[len(ws_) // 2], ws_[0], sum(1 for x in nl_ if x == 1),
              coherence([mem[k] for k in alive if k not in (p, q)] + [mem[p] + mem[q]])))
        for k in sorted([k for k in alive if k not in (p, q)] + [r_], key=lambda k: -w[k]): print('   [%s] %d particles  %s' % (' '.join(mem[k]) if len(mem[k]) <= 12 else '%d leaves' % len(mem[k]), w[k], describe(mem[k])))
    steps.append((len(steps) + 1, p, q, r_, z, d, wp, wq)); alive = [k for k in alive if k not in (p, q)] + [r_]
    print('merge %2d: d2 %.4f z %7.1f  [%s] (%d) + [%s] (%d)%s' % (len(steps), d, z, ' '.join(mem[p]) if len(mem[p]) <= 12 else '%d leaves' % len(mem[p]), wp, ' '.join(mem[q]) if len(mem[q]) <= 12 else '%d leaves' % len(mem[q]), wq,
          ('\n      %s\n    + %s' % (describe(mem[p]), describe(mem[q]))) if (conf is not None or arm) else ''))
if cut is None: cut = [mem[alive[0]]]
nw_ = {i: names[i] for i in range(N)}; ht = {i: 0. for i in range(N)}   # the tree in Newick, branch lengths = merge height (d2 or z) differences
for st_, p_, q_, r_, z_, d_, _, _ in steps:
    h_ = d_ if a.by == 'd2' else z_; ht[r_] = h_; nw_[r_] = '(%s:%.4f,%s:%.4f)' % (nw_[p_], max(h_ - ht[p_], 0), nw_[q_], max(h_ - ht[q_], 0))
open(a.out + '/tree.nwk', 'w').write(nw_[steps[-1][3]] + ';\n' if steps else '')
if a.by == 'd2': cut = []   # by d2 the tree is the result; cut_groups only for --by z
with open(a.out + '/merge_tree.txt', 'w') as fo:
    fo.write('# step node_a node_b new_node z d2 count_a count_b (merged by %s; nodes < %d are the leaves %s)\n' % (a.by, N, ' '.join(names)))
    for s in steps: fo.write('%d %d %d %d %.2f %.5f %d %d\n' % s)
if cut: print('== cut at z %.1f: %d distinguishable states' % (a.zthr, len(cut)))
with open(a.out + '/cut_groups.txt', 'w') as fo:
    for g in sorted(cut, key=lambda g: -sum(cnt[names.index(m)] for m in g)):
        nn = sum(cnt[names.index(m)] for m in g); fo.write('%s\n' % ' '.join(g)); print('   [%s] %d particles  %s' % (' '.join(g), nn, describe(g)))

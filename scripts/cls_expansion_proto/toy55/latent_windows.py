#!/usr/bin/env python3
"""Coarse-to-fine latent windows on the RNase P toy: a continuous particle coordinate -> equal-population windows ->
window reconstructions (SIMPLE, outside) -> soft cross-half E-step of every particle against the window maps with a
smoothness prior along the coordinate -> posterior-mean coordinate for the next round. Everything scored against the
truth conformers. Subcommands:
  assign  COORD GT --nptcl N --nw 8 --out states.txt --info windows.txt --states .. --axis ..   [--col 1]
  estep   STK TABLE GT --mapdir DIR --nw 8 --lp 14 --out coord.txt --states .. --axis .. [--prev states.txt]
          [--prior_sigma 1.5] [--refs f1,f2,.. (oracle: truth maps in --states order, no halves)] [--nsub 0] [--info windows.txt]
  maps    MAPDIR --nw 8 --info windows.txt --truth f1,f2,.. --states .. --cons toy_consensus.mrc
  merge   MASTER STATES --tables w1:file1,w2:file2,.. --out merged.txt   (per-window refined poses back into one table)
TABLE = print_project_field ptcl3D (e1 e2 e3 x y dfx dfy angast eo state), pind = line number."""
import sys, os, argparse, time
import numpy as np
os.environ.setdefault('TOY_SMPD', '1.464')
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/pvbag/subclass_pool'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from view_tangent_toy import euler2m, ctf_eval, Band, fft2c, SMPD
from scipy.ndimage import map_coordinates
from scipy.stats import spearmanr
from scipy.optimize import linear_sum_assignment
from model_class import read_mrc, Shells, fftc, BANDS


def read_gt(fn):
    gt = {}
    for l in open(fn):
        if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = (int(w[1]), int(w[2]))
    return gt


def read_table(fn, keys):
    rows = []
    for l in open(fn):
        if not l.strip() or l.startswith('#'): continue
        kv = dict(t.split('=', 1) for t in l.split() if '=' in t); rows.append([float(kv.get(k, -1 if k == 'eo' else 0)) for k in keys])
    return np.array(rows)


def read_coord(fn, col=1):
    out = {}
    for l in open(fn):
        if l.strip() and not l.startswith('#'): w = l.split(); out[int(w[0])] = float(w[col])
    return out


def hungarian_acc(pred, truth):
    pu, pi = np.unique(pred, return_inverse=True); tu, ti = np.unique(truth, return_inverse=True)
    C = np.zeros((len(pu), len(tu))); np.add.at(C, (pi, ti), 1); r, c = linear_sum_assignment(-C); return C[r, c].sum() / len(pred)


def readout_windows(win, pinds, gt, states, axis, nw, tag, info=None, extra=''):
    conf = np.array([gt.get(int(p), (0, 0))[0] for p in pinds]); closed = np.array([gt.get(int(p), (0, 0))[1] for p in pinds]) == 1; good = conf > 0
    ax = np.array([axis.get(c, 0.) for c in conf]); lines = []; pur = 0.; best_closed = (0, 0, 0)
    for w in range(1, nw + 1):
        m = (win == w); mg = m & good; n = m.sum()
        if mg.sum() == 0: lines.append('%d %d 0 0.00 0.00 0.00 0.00' % (w, n)); continue
        h = np.bincount(conf[mg]); dom = h.argmax(); share = h[dom] / mg.sum(); pur += h[dom]
        cs = closed[mg].mean(); lines.append('%d %d %d %.3f %.3f %.3f %.3f' % (w, n, dom, share, cs, ax[mg].mean(), 1 - mg.sum() / max(n, 1)))
        nc = np.sum(closed & m)
        if nc > best_closed[0]: best_closed = (nc, w, nc / max(n, 1))
    rho = spearmanr(win[good], ax[good])[0] if good.sum() > 2 else 0.
    print('%s windows (%d): Spearman(window, truth axis) %+.3f | weighted purity %.3f | Hungarian accuracy %.3f | closed recall %.2f at purity %.2f (window %d)%s' %
          (tag, nw, rho, pur / max(good.sum(), 1), hungarian_acc(win[good], conf[good]), best_closed[0] / max(closed.sum(), 1), best_closed[2], best_closed[1], extra))
    for l in lines:
        w, n, dom, share, cs, am, js = l.split(); print('%s   window %s: n %6s, dominant conformer %2s (share %s), closed %s, truth axis mean %s, junk %s' % (tag, w, n, dom, share, cs, am, js))
    if info:
        with open(info, 'w') as f: f.write('# window n dominant share closed_share axis_mean junk_share\n' + '\n'.join(lines) + '\n')


def cmd_assign(a):
    coord = read_coord(a.coord, a.col); gt = read_gt(a.gt); axis = parse_axis(a)
    pinds = np.array(sorted(coord)); c = np.array([coord[p] for p in pinds])
    # rank-based equal-population windows (random tie-breaking: ties cannot empty a window)
    rng = np.random.default_rng(0); order = np.lexsort((rng.random(len(c)), c)); rank = np.empty(len(c), int); rank[order] = np.arange(len(c)); win = (rank * a.nw) // len(c) + 1
    edges = np.quantile(c, np.linspace(0, 1, a.nw + 1))
    states = np.zeros(a.nptcl, int); states[pinds - 1] = win; np.savetxt(a.out, states, fmt='%d')
    readout_windows(win, pinds, gt, a.states, axis, a.nw, a.tag, a.info, extra=' | edges %s' % np.array2string(edges, precision=2))


def parse_axis(a):
    states = [int(s) for s in a.states.split(',')]; return dict(zip(states, [float(v) for v in a.axis.split(',')]))


def slice_batch(VF_re, VF_im, Rs, kx, ky, conv):
    c = VF_re.shape[0] // 2; hk = np.stack([kx, ky, np.zeros_like(kx)], 1); M = Rs if conv == 0 else np.transpose(Rs, (0, 2, 1))
    loc = np.einsum('kj,mji->mki', hk, M); coords = [loc[..., 2].ravel() + c, loc[..., 1].ravel() + c, loc[..., 0].ravel() + c]
    re = map_coordinates(VF_re, coords, order=1, mode='nearest'); im = map_coordinates(VF_im, coords, order=1, mode='nearest')
    return (re + 1j * im).reshape(len(Rs), -1)


class Data:
    """band-limited Fourier data of a particle chunk at the particle poses; stk = one memmap, or a list of memmaps with
    per-particle (stack id, index) arrays sid/iid; ctfabs = |CTF| model (phase-flipped stacks)"""
    def __init__(self, stk, tab, box, Bd, ix, iy, sign, sid=None, iid=None, ctfabs=False):
        self.stk, self.tab, self.box, self.Bd, self.ix, self.iy, self.sign, self.sid, self.iid, self.ctfabs = stk, tab, box, Bd, ix, iy, sign, sid, iid, ctfabs

    def image(self, i):
        return np.array(self.stk[self.sid[i]][self.iid[i]], float) if self.sid is not None else np.array(self.stk[i], float)

    def chunk(self, idx):
        F = np.empty((len(idx), self.Bd.n), complex); C = np.empty((len(idx), self.Bd.n)); R = np.empty((len(idx), 3, 3))
        for t, i in enumerate(idx):
            im = self.image(i); Fi = fft2c(im); F[t] = Fi[self.iy, self.ix] * self.sign
            C[t] = ctf_eval(self.Bd.kx, self.Bd.ky, self.tab[i, 5], self.tab[i, 6], self.tab[i, 7], self.box); R[t] = euler2m([-self.tab[i, 0], self.tab[i, 1], -self.tab[i, 2]])
        if self.ctfabs: C = np.abs(C)
        return F, C, R


def shift_phase(Bd, sx, sy, box, ss, fs):
    return np.exp(-2j * np.pi * (Bd.kx[None, :] * (ss * fs * sx)[:, None] + Bd.ky[None, :] * (ss * fs * sy)[:, None]) / box)


def calibrate(stk, tab, box, Bd, ix, iy, VF, cal_idx, sid=None, iid=None, ctfabs=False):
    """projection convention, shift sign/scale and contrast sign on a subset against the consensus FT VF (re, im)"""
    D = Data(stk, tab, box, Bd, ix, iy, 1., sid=sid, iid=iid, ctfabs=ctfabs); F, C, R = D.chunk(cal_idx); best = None
    for conv in (0, 1):
        P0 = slice_batch(VF[0], VF[1], R, Bd.kx, Bd.ky, conv)
        for ss in (1., -1.):
            for fs in (1.0, 0.5):
                P = C * P0 * shift_phase(Bd, tab[cal_idx, 3], tab[cal_idx, 4], box, ss, fs)
                cc = np.real(np.sum(Bd.wq * np.conj(F) * P, 1)) / np.sqrt(np.sum(Bd.wq * np.abs(F) ** 2, 1) * np.sum(Bd.wq * np.abs(P) ** 2, 1)); m = float(np.mean(cc))
                if best is None or abs(m) > abs(best[0]): best = (m, conv, ss, fs)
    m, conv, ss, fs = best; sign = 1. if m > 0 else -1.
    print('   calibration: conv %d, shift sign %+.0f, scale %.1f, contrast sign %+.0f (mean ncc %.3f on %d particles)' % (conv, ss, fs, sign, abs(m), len(cal_idx)))
    return conv, ss, fs, sign


def load_refs(a, S, msk):
    """returns refs[half][w] = (re, im) of the masked map FT; half 0 = even map, 1 = odd map (used cross-wise)"""
    if a.refs:
        fns = a.refs.split(','); M = [read_mrc(f) * msk for f in fns]; FT = [fftc(v) for v in M]
        R = [(np.ascontiguousarray(f.real, np.float32), np.ascontiguousarray(f.imag, np.float32)) for f in FT]
        return [R, R], len(fns), np.mean(M, 0)
    refs = [[], []]; full = []
    for w in range(1, a.nw + 1):
        for h, suf in enumerate(('_even', '_odd')):
            fn = '%s/recvol_state%02d%s.mrc' % (a.mapdir, w, suf)
            if not os.path.exists(fn): fn = '%s/recvol_state%02d.mrc' % (a.mapdir, w)
            f = fftc(read_mrc(fn) * msk); refs[h].append((np.ascontiguousarray(f.real, np.float32), np.ascontiguousarray(f.imag, np.float32)))
        full.append(read_mrc('%s/recvol_state%02d.mrc' % (a.mapdir, w)) * msk)
    return refs, a.nw, np.mean(full, 0)


def cmd_estep(a):
    t0 = time.time(); gt = read_gt(a.gt); axis = parse_axis(a); states = [int(s) for s in a.states.split(',')]
    tab = read_table(a.table, ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast', 'eo', 'state')); N = len(tab)
    h = np.fromfile(a.stk, np.int32, 256); box = int(h[0]); stk = np.memmap(a.stk, np.float32, 'r', offset=1024 + int(h[23]), shape=(int(h[2]), box, box))
    eo = tab[:, 8].astype(int); noeo = eo < 0; eo[noeo] = np.arange(N)[noeo] % 2
    prev = np.loadtxt(a.prev).astype(int) if a.prev else np.zeros(N, int)
    sel = np.where((tab[:, 9] > 0) & ((prev > 0) if a.prev else True))[0]
    if a.nsub and len(sel) > a.nsub: sel = np.sort(np.random.default_rng(0).choice(sel, a.nsub, replace=False))
    S = Shells(box, SMPD); msk = S.soft_mask(a.mskdiam / 2 / SMPD); refs, K, cons = load_refs(a, S, msk)
    kmax = int(np.floor(box * SMPD / a.lp)); Bd = Band(box, kmax); c = box // 2; ix = (Bd.kx + c).astype(int); iy = (Bd.ky + c).astype(int)
    print('%s E-step: %d particles, %d references%s, band to %.1f A (%d shells, %d coefficients)' % (a.tag, len(sel), K, '' if a.refs else ' (cross-half)', a.lp, kmax, Bd.n), flush=True)
    rng = np.random.default_rng(1); conf_all = np.array([gt.get(i + 1, (0, 0))[0] for i in range(N)])
    cal = rng.choice(sel[conf_all[sel] > 0], min(300, len(sel)), replace=False); VC = fftc(cons); VC = (np.ascontiguousarray(VC.real, np.float32), np.ascontiguousarray(VC.imag, np.float32))
    conv, ss, fs, sign = calibrate(stk, tab, box, Bd, ix, iy, VC, cal)
    D = Data(stk, tab, box, Bd, ix, iy, sign)
    # noise per shell from the residual against the consensus
    sig_idx = rng.choice(sel, min(2000, len(sel)), replace=False); F, C, R = D.chunk(sig_idx); P = C * slice_batch(VC[0], VC[1], R, Bd.kx, Bd.ky, conv) * shift_phase(Bd, tab[sig_idx, 3], tab[sig_idx, 4], box, ss, fs)
    # per-shell amplitude scale of the references to the images (external maps are in arbitrary units; without it the
    # likelihood ranks references by their own power), least squares against the consensus on the subset, smoothed, floored
    shells = [Bd.shell == s for s in range(kmax + 1)]
    num = np.array([np.sum(Bd.wq[m] * np.real(np.conj(P[:, m]) * F[:, m])) if m.any() else 0. for m in shells]); den = np.array([np.sum(Bd.wq[m] * np.abs(P[:, m]) ** 2) if m.any() else 0. for m in shells])
    sc = np.where(den > 0, num / np.maximum(den, 1e-30), 0.); sc = np.convolve(np.pad(sc, 1, mode='edge'), np.ones(3) / 3., mode='valid'); sc = np.maximum(sc, 0.05 * max(sc[1:5].mean(), 1e-30))
    SC = sc[Bd.shell]; P = P * SC
    res2 = np.abs(F - P) ** 2; sig2 = np.array([res2[:, m].mean() if m.any() else 1. for m in shells]); sig2 = np.maximum(sig2, 1e-12)[Bd.shell]
    snr = np.array([np.mean(np.abs(P[:, Bd.shell == s]) ** 2) / sig2[Bd.shell == s][0] for s in range(1, kmax + 1)])
    print('   reference scale per shell (relative to shell 1): %s' % np.array2string(sc[1:: max(1, kmax // 8)] / sc[1], precision=2), flush=True)
    print('   model power / noise per shell (signal fraction): %s' % np.array2string(snr[:: max(1, kmax // 8)], precision=3), flush=True)
    ll = np.zeros((len(sel), K)); CH = 1000
    for s0 in range(0, len(sel), CH):
        idx = sel[s0:s0 + CH]; F, C, R = D.chunk(idx); SH = shift_phase(Bd, tab[idx, 3], tab[idx, 4], box, ss, fs); CS = C * SH
        for half in (0, 1):
            m = eo[idx] == half
            if not m.any(): continue
            for w in range(K):
                re, im = refs[1 - half][w]; P = CS[m] * slice_batch(re, im, R[m], Bd.kx, Bd.ky, conv) * SC
                ll[s0:s0 + CH][m, w] = -np.sum(Bd.wq * np.abs(F[m] - P) ** 2 / sig2, 1)
        if s0 % (10 * CH) == 0: print('   %d / %d particles (%.0f s)' % (s0 + len(idx), len(sel), time.time() - t0), flush=True)
    np.save(a.out + '.ll.npy', ll.astype(np.float32))
    conf = conf_all[sel]; closed = np.array([gt.get(i + 1, (0, 0))[1] for i in sel]) == 1; good = conf > 0; ax = np.array([axis.get(c, 0.) for c in conf])
    wi = np.arange(1, K + 1)
    def posterior(logprior):
        L = ll + logprior; L -= L.max(1, keepdims=True); p = np.exp(L); return p / p.sum(1, keepdims=True)
    outs = {'free': posterior(np.zeros((len(sel), K)))}
    if a.prev and not a.refs: outs['smooth'] = posterior(-0.5 * ((wi[None, :] - prev[sel][:, None]) / a.prior_sigma) ** 2)
    for name, post in outs.items():
        am = post.argmax(1) + 1; mean = post @ wi; mx = post.max(1); ent = -np.sum(post * np.log(post + 1e-12), 1) / np.log(K)
        if a.refs:    # oracle: references are the truth conformers in --states order
            pred = np.array([states[i] for i in post.argmax(1)]); acc = np.mean(pred[good] == conf[good]); axm = post @ np.array([axis[s] for s in states])
            print('%s oracle (%s prior): per-particle accuracy %.3f over %d conformer particles (chance %.3f); mean max posterior %.2f; Spearman(posterior-mean axis, truth) %+.3f; closed recall %.2f' %
                  (a.tag, name, acc, good.sum(), np.bincount(conf[good]).max() / good.sum(), mx[good].mean(), spearmanr(axm[good], ax[good])[0], np.mean(pred[good & closed] == conf[good & closed]) if (good & closed).any() else 0))
            for s in states: m = good & (conf == s); print('%s   conformer %2d: recall %.3f (n %d)' % (a.tag, s, np.mean(pred[m] == s) if m.any() else 0, m.sum()))
        else:
            moved = np.mean(am != prev[sel]) if a.prev else 0.
            print('%s E-step (%s prior): Spearman(posterior mean, truth axis) %+.3f | mean max posterior %.2f, entropy %.2f | moved window %.2f' % (a.tag, name, spearmanr(mean[good], ax[good])[0], mx[good].mean(), ent[good].mean(), moved))
            readout_windows(am, sel + 1, gt, a.states, axis, K, a.tag + ' ' + name)
    post = outs['free']; am = post.argmax(1) + 1; mean = post @ wi      # the coordinate carried forward is the pure-EM (flat prior) posterior; the smooth one is a readout
    with open(a.out, 'w') as f:
        f.write('# pind coord_mean argmax maxpost entropy\n')
        for p, m_, a_, x_, e_ in zip(sel + 1, mean, am, post.max(1), -np.sum(post * np.log(post + 1e-12), 1) / np.log(K)): f.write('%d %.5f %d %.4f %.4f\n' % (p, m_, a_, x_, e_))
    np.savez(a.out + '.post.npz', pinds=sel + 1, post=post.astype(np.float32), calib=np.array([conv, ss, fs, sign]))
    print('%s posteriors -> %s (%.0f s)' % (a.tag, a.out, time.time() - t0))


def cmd_maps(a):
    axis = parse_axis(a); states = [int(s) for s in a.states.split(',')]; T = {s: read_mrc(f) for s, f in zip(states, a.truth.split(','))}
    cons = read_mrc(a.cons); B = cons.shape[0]; S = Shells(B, SMPD); msk = S.soft_mask(a.mskdiam / 2 / SMPD)
    info = np.loadtxt(a.info, comments='#', ndmin=2); n = info[:, 1]; dom = info[:, 2].astype(int)
    M = [read_mrc('%s/recvol_state%02d.mrc' % (a.mapdir, w)) * msk for w in range(1, a.nw + 1)]; w_ = n / n.sum(); Mm = sum(wi * m for wi, m in zip(w_, M))
    FT = {s: fftc((T[s] - cons) * msk) for s in states}; sh = S.sh
    def fsc_res(x, y):
        a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
        f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
        while s < B // 2 and f[s] > 0.143: s += 1
        return B * SMPD / s
    def bc(x, y, lo, hi):
        m = S.band(lo, hi); return float(np.real(np.sum(np.conj(x[m]) * y[m])) / max(np.sqrt(np.sum(np.abs(x[m]) ** 2) * np.sum(np.abs(y[m]) ** 2)), 1e-30))
    right = 0
    for w in range(1, a.nw + 1):
        D = fftc(M[w - 1] - Mm); d = dom[w - 1]
        if d == 0: print('%s   window %d: no conformer particles' % (a.tag, w)); continue
        cc = [bc(D, FT[d], lo, hi) for lo, hi in BANDS]; best = max(states, key=lambda s: bc(D, FT[s], 20, 8)); right += best == d
        try: res = fsc_res(read_mrc('%s/recvol_state%02d_even.mrc' % (a.mapdir, w)) * msk, read_mrc('%s/recvol_state%02d_odd.mrc' % (a.mapdir, w)) * msk)
        except Exception: res = 0.
        print('%s   window %d (n %d, dominant %d at %.2f): FSC 0.143 %.1f A | corr(map - mean, truth %d - cons) per band %s | best conformer at 20-8 A: %d %s' %
              (a.tag, w, n[w - 1], d, info[w - 1, 3], res, d, ' '.join('%.2f' % c for c in cc), best, 'OK' if best == d else 'MISMATCH'))
    print('%s maps: %d / %d windows match their dominant conformer at 20-8 A' % (a.tag, right, a.nw))


def cmd_merge(a):
    master = [l.rstrip('\n') for l in open(a.master) if l.strip() and not l.startswith('#')]; states = np.loadtxt(a.states).astype(int); N = len(master)
    keys = ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'angast'); out = [None] * N; nrep = 0
    def fmt(line, st):
        kv = dict(t.split('=', 1) for t in line.split() if '=' in t); return ' '.join('%s=%s' % (k, kv.get(k, '0')) for k in keys) + ' state=%d' % st
    for i in range(N): out[i] = fmt(master[i], states[i])
    for item in a.tables.split(','):
        w, fn = item.split(':'); w = int(w); rows = [l.rstrip('\n') for l in open(fn) if l.strip() and not l.startswith('#')]; assert len(rows) == N, (fn, len(rows), N)
        for i in np.where(states == w)[0]: out[i] = fmt(rows[i], w); nrep += 1
    open(a.out, 'w').write('\n'.join(out) + '\n'); print('merge: %d particles, %d rows replaced from %d window tables -> %s' % (N, nrep, len(a.tables.split(',')), a.out))


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest='cmd')
    p = sub.add_parser('assign'); p.add_argument('coord'); p.add_argument('gt'); p.add_argument('--nptcl', type=int, required=True); p.add_argument('--nw', type=int, default=8); p.add_argument('--out', required=True)
    p.add_argument('--info', default=''); p.add_argument('--states', required=True); p.add_argument('--axis', required=True); p.add_argument('--col', type=int, default=1); p.add_argument('--tag', default='assign')
    p = sub.add_parser('estep'); p.add_argument('stk'); p.add_argument('table'); p.add_argument('gt'); p.add_argument('--mapdir', default=''); p.add_argument('--refs', default=''); p.add_argument('--nw', type=int, default=8)
    p.add_argument('--lp', type=float, default=14.); p.add_argument('--out', required=True); p.add_argument('--states', required=True); p.add_argument('--axis', required=True); p.add_argument('--prev', default='')
    p.add_argument('--prior_sigma', type=float, default=1.5); p.add_argument('--nsub', type=int, default=0); p.add_argument('--mskdiam', type=float, default=200.); p.add_argument('--tag', default='estep')
    p = sub.add_parser('maps'); p.add_argument('mapdir'); p.add_argument('--nw', type=int, default=8); p.add_argument('--info', required=True); p.add_argument('--truth', required=True); p.add_argument('--states', required=True)
    p.add_argument('--axis', required=True); p.add_argument('--cons', required=True); p.add_argument('--mskdiam', type=float, default=200.); p.add_argument('--tag', default='maps')
    p = sub.add_parser('merge'); p.add_argument('master'); p.add_argument('states'); p.add_argument('--tables', required=True); p.add_argument('--out', required=True)
    a = ap.parse_args(); {'assign': cmd_assign, 'estep': cmd_estep, 'maps': cmd_maps, 'merge': cmd_merge}[a.cmd](a)


if __name__ == '__main__':
    main()

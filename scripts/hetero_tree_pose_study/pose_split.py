#!/usr/bin/env python3
"""pose_split.py: pose-change statistics of every candidate split of saved hetero_tree runs, read only.
For every node with n.simple (the node's input poses) and stageB/b.simple (its refined poses and states): per particle the angle
between the input and the stage-B rotation (symmetry-reduced for cN) and the change of the in-plane shift; per candidate pair of
children (the delivered split groups of split/split_children.txt where the node split, the stage-B states where it did not) the
flip fractions (> 30 deg) of the two children, their difference, the KS distance between the two children's angle distributions and
between their shift-change distributions, each with a z against NSHUF label shuffles inside the split test's nuisance strata
(40 Fibonacci view bins of the input pose x defocus tertile), plus the children's view-distribution distance (TV over the 40 bins)
as the control statistic. One TSV row per pair and one per node.
usage: pose_split.py OUT_PREFIX NAME=TREE_DIR[:PGRP] ...   (TREE_DIR holds 1_hetero_tree or the node_* folders; PGRP c1 default)
env NPROC (processes), NSHUF (shuffles, 200), BUILD (an install whose simple_exec reads the projects)"""
import glob, os, subprocess, sys
import numpy as np
from multiprocessing import Pool

KEYS = ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'state')
NVIEW, NSHUF, FLIP = 40, int(os.environ.get('NSHUF', 200)), 30.


def dump(proj):
    p = subprocess.run(['simple_exec', 'prg=print_project_field', 'projfile=' + proj, 'oritype=ptcl3D', 'mkdir=no'],
                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    out = []
    for ln in p.stdout.splitlines():
        if 'e1=' not in ln:
            continue
        d = {}
        for t in ln.split():
            k, _, v = t.partition('=')
            if k in KEYS:
                d[k] = float(v)
        out.append([d.get(k, 0.) for k in KEYS])
    return np.array(out).reshape(-1, len(KEYS))


def rz(t):
    R = np.zeros((len(t), 3, 3)); R[:, 0, 0] = np.cos(t); R[:, 0, 1] = -np.sin(t); R[:, 1, 0] = np.sin(t); R[:, 1, 1] = np.cos(t); R[:, 2, 2] = 1
    return R


def ry(t):
    R = np.zeros((len(t), 3, 3)); R[:, 0, 0] = np.cos(t); R[:, 0, 2] = -np.sin(t); R[:, 2, 0] = np.sin(t); R[:, 2, 2] = np.cos(t); R[:, 1, 1] = 1
    return R


def mats(e):  # SIMPLE's euler2m: R = Rz(e3) Ry(e2) Rz(e1)
    a, b, c = np.radians(e[:, 0]), np.radians(e[:, 1]), np.radians(e[:, 2])
    return rz(c) @ ry(b) @ rz(a)


def ang(e0, e1, nsym):  # angle between R0 and R1 S_k, minimised over the cN group (equivalent poses of a symmetric volume)
    R0, R1 = mats(e0), mats(e1)
    best = None
    for k in range(nsym):
        S = rz(np.full(1, 2 * np.pi * k / nsym))[0]
        E = np.einsum('nji,njk->nik', R0, R1 @ S)
        a = np.degrees(np.arccos(np.clip((np.trace(E, axis1=1, axis2=2) - 1) / 2, -1, 1)))
        best = a if best is None else np.minimum(best, a)
    return best


def fib_dirs():
    i = np.arange(NVIEW); t = i + 0.5; z = 1 - 2 * t / NVIEW; r = np.sqrt(np.maximum(1 - z * z, 0)); g = np.pi * (3 - np.sqrt(5)) * t
    return np.stack([r * np.cos(g), r * np.sin(g), z], 1)


def view_bins(e):
    e1, e2 = np.radians(e[:, 0]), np.radians(e[:, 1])
    v = np.stack([np.sin(e2) * np.cos(e1), np.sin(e2) * np.sin(e1), np.cos(e2)], 1)
    return np.argmax(v @ fib_dirs().T, 1)


def strata(e, df):
    q = np.quantile(df, [1 / 3, 2 / 3]) if len(df) > 2 else np.array([np.inf, np.inf])
    return view_bins(e) * 3 + (df > q[0]).astype(int) + (df > q[1]).astype(int)


def ks(a, b):
    if len(a) == 0 or len(b) == 0:
        return 0.
    x = np.sort(np.concatenate([a, b])); ca = np.searchsorted(np.sort(a), x, 'right') / len(a); cb = np.searchsorted(np.sort(b), x, 'right') / len(b)
    return float(np.max(np.abs(ca - cb)))


def tv_views(va, vb):
    ha = np.bincount(va, minlength=NVIEW) / max(len(va), 1); hb = np.bincount(vb, minlength=NVIEW) / max(len(vb), 1)
    return 0.5 * float(np.abs(ha - hb).sum())


def pair_stats(sel_a, sel_b, dang, dsh, st, rng):
    fa, fb = float(np.mean(dang[sel_a] > FLIP)), float(np.mean(dang[sel_b] > FLIP))
    obs = np.array([abs(fa - fb), ks(dang[sel_a], dang[sel_b]), ks(dsh[sel_a], dsh[sel_b])])
    idx = np.concatenate([sel_a, sel_b]); lab = np.concatenate([np.zeros(len(sel_a), bool), np.ones(len(sel_b), bool)])
    s = st[idx]; order = np.argsort(s, kind='stable'); s_sorted = s[order]; bounds = np.flatnonzero(np.diff(s_sorted)) + 1
    groups = np.split(order, bounds)
    nul = np.zeros((NSHUF, 3))
    for r in range(NSHUF):
        pl = lab.copy()
        for g in groups:
            if len(g) > 1:
                pl[g] = lab[g][rng.permutation(len(g))]
        a, b = idx[~pl], idx[pl]
        nul[r] = [abs(np.mean(dang[a] > FLIP) - np.mean(dang[b] > FLIP)), ks(dang[a], dang[b]), ks(dsh[a], dsh[b])]
    mu, sd = nul.mean(0), nul.std(0, ddof=1)
    z = (obs - mu) / np.maximum(sd, 1e-9)
    return fa, fb, obs, z, mu, sd


def node(job):
    name, nd, nsym = job
    nid = os.path.basename(nd)[5:]
    rng = np.random.default_rng(1 + len(nid))
    n0, nb = dump(nd + '/n.simple'), dump(nd + '/stageB/b.simple')
    rows = np.loadtxt(nd + '/rows.txt', dtype=int, ndmin=1)
    if len(n0) == 0 or len(nb) == 0:
        return [], ['\t'.join([name, nid, str(len(nid) - 1), '0', '-', '-', '-', '-', 'no dump'])]
    if len(n0) != len(rows):  # master-sized projects: the node's rows index them
        n0, nb = n0[rows - 1], nb[rows - 1]
    k = min(len(n0), len(nb)); n0, nb = n0[:k], nb[:k]
    act = nb[:, 7] > 0
    dang = np.full(k, np.nan); dsh = np.full(k, np.nan)
    dang[act] = ang(n0[act, :3], nb[act, :3], nsym)
    dsh[act] = np.hypot(nb[act, 3] - n0[act, 3], nb[act, 4] - n0[act, 4])
    st = strata(n0[:, :3], 0.5 * (n0[:, 5] + n0[:, 6]))
    vb = view_bins(n0[:, :3])
    # stA0 flips for the node line, when the step's project survived pruning
    f_a0 = '-'
    if os.path.exists(nd + '/stA0/r.simple'):
        a0 = dump(nd + '/stA0/r.simple')
        if len(a0) != len(rows) and len(a0) >= rows.max():
            a0 = a0[rows - 1]
        m = min(len(a0), k); am = a0[:m, 7] > 0
        if am.sum() > 0:
            f_a0 = '%.4f' % np.mean(ang(n0[:m][am, :3], a0[:m][am, :3], nsym) > FLIP)
    # candidate pairs: the delivered groups where the node split, else the stage-B states
    kind, lab = 'refused', nb[:, 7].astype(int)
    sc = nd + '/split/split_children.txt'
    if os.path.exists(sc):
        ch = np.loadtxt(sc, dtype=int, ndmin=1)
        if len(ch) == k and len(np.unique(ch[ch > 0])) >= 2:
            kind, lab = 'split', ch
    glab = np.unique(lab[(lab > 0) & act])
    pairs = []
    for i in range(len(glab)):
        for j in range(i + 1, len(glab)):
            a, b = np.flatnonzero((lab == glab[i]) & act), np.flatnonzero((lab == glab[j]) & act)
            if len(a) < 50 or len(b) < 50:
                continue
            fa, fb, obs, z, mu, sd = pair_stats(a, b, dang, dsh, st, rng)
            pairs.append('\t'.join([name, nid, kind, '%d|%d' % (glab[i], glab[j]), str(len(a)), str(len(b)),
                                    '%.4f' % fa, '%.4f' % fb, '%.4f' % obs[0], '%.2f' % z[0], '%.4f' % mu[0], '%.4f' % sd[0],
                                    '%.4f' % obs[1], '%.2f' % z[1], '%.4f' % obs[2], '%.2f' % z[2],
                                    '%.2f' % np.median(dang[a]), '%.2f' % np.median(dang[b]), '%.3f' % tv_views(vb[a], vb[b])]))
    nl = '\t'.join([name, nid, str(len(nid) - 1), str(int(act.sum())), f_a0, '%.4f' % np.mean(dang[act] > FLIP),
                    '%.4f' % np.mean(dang[act] > 10), '%.2f' % np.median(dang[act]), kind + ' %d groups' % len(glab)])
    return pairs, [nl]


def selftest():
    e = np.array([[30., 60., 20.], [100., 120., 300.]])
    assert np.allclose(ang(e, e, 1), 0, atol=1e-4)
    assert np.allclose(ang(e, e + [180, 0, 0], 2), 0, atol=1e-3) and np.allclose(ang(e, e + [180, 0, 0], 1), 180, atol=1e-3)
    rng = np.random.default_rng(0); d = np.abs(rng.normal(0, 5, 2000)); d[:300] = 170; dsh = rng.normal(0, 1, 2000); st = rng.integers(0, 120, 2000)
    fa, fb, obs, z, mu, sd = pair_stats(np.arange(0, 1000), np.arange(1000, 2000), d, dsh, st, rng)
    assert z[0] > 5 and z[1] > 5, (obs, z)
    pm = rng.permutation(2000); fa, fb, obs, z, mu, sd = pair_stats(pm[:1000], pm[1000:], d, dsh, st, rng)
    assert abs(z[0]) < 4 and abs(z[1]) < 4, (obs, z)
    print('selftest ok: artefact pair z', z)


if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        selftest(); sys.exit(0)
    out = sys.argv[1]
    jobs = []
    for arg in sys.argv[2:]:
        name, d = arg.split('=', 1)
        pg = 'c1'
        if ':' in d:
            d, pg = d.rsplit(':', 1)
        nsym = int(pg[1:]) if pg[0] == 'c' and pg[1:].isdigit() else 1
        T = d + '/1_hetero_tree' if os.path.isdir(d + '/1_hetero_tree') else d
        for nd in sorted(glob.glob(T + '/node_*')):
            if os.path.exists(nd + '/n.simple') and os.path.exists(nd + '/stageB/b.simple'):
                jobs.append((name, nd, nsym))
    with Pool(int(os.environ.get('NPROC', 8))) as p:
        res = p.map(node, jobs, chunksize=1)
    with open(out + '_pairs.tsv', 'w') as f:
        f.write('tree\tnode\tkind\tpair\tn_a\tn_b\tflip_a\tflip_b\tdflip\tz_dflip\tnull_dflip_mu\tnull_dflip_sd\tks_ang\tz_ks_ang\tks_shift\tz_ks_shift\tmed_a\tmed_b\tview_tv\n')
        for pr, _ in res:
            f.writelines(r + '\n' for r in pr)
    with open(out + '_nodes.tsv', 'w') as f:
        f.write('tree\tnode\tdepth\tn_refined\tflip30_stA0\tflip30_stB\tflip10_stB\tmed_stB\tcandidates\n')
        for _, nl in res:
            f.writelines(r + '\n' for r in nl)
    print('wrote', out, len(jobs), 'nodes,', sum(len(pr) for pr, _ in res), 'pairs')

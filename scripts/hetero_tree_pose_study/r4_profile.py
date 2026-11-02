#!/usr/bin/env python3
"""r4_profile.py TREE_DIR OUT.txt: the 692k root's four groups and r4's two groups compared per particle, read only.
Per group: size, the pose change the node's refinement made (fraction > 30, > 90, > 150 degrees, median), the stage-B corr
(quality proxy) quantiles, mean defocus, the odd/even balance, the view distribution's TV distance to the rest of the node,
and the same pose-change fractions split by corr quartile (do the turned particles come from the low-quality end?)."""
import os, subprocess, sys
import numpy as np

KEYS = ('e1', 'e2', 'e3', 'x', 'y', 'dfx', 'dfy', 'state', 'corr', 'eo')


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


def mats(e):
    a, b, c = np.radians(e[:, 0]), np.radians(e[:, 1]), np.radians(e[:, 2])
    return rz(c) @ ry(b) @ rz(a)


def ang(e0, e1):
    E = np.einsum('nji,njk->nik', mats(e0), mats(e1))
    return np.degrees(np.arccos(np.clip((np.trace(E, axis1=1, axis2=2) - 1) / 2, -1, 1)))


def views(e, nv=40):
    i = np.arange(nv); t = i + 0.5; z = 1 - 2 * t / nv; r = np.sqrt(np.maximum(1 - z * z, 0)); g = np.pi * (3 - np.sqrt(5)) * t
    F = np.stack([r * np.cos(g), r * np.sin(g), z], 1)
    e1, e2 = np.radians(e[:, 0]), np.radians(e[:, 1])
    v = np.stack([np.sin(e2) * np.cos(e1), np.sin(e2) * np.sin(e1), np.cos(e2)], 1)
    return np.argmax(v @ F.T, 1)


def profile(nd, out, title):
    n0, nb = dump(nd + '/n.simple'), dump(nd + '/stageB/b.simple')
    k = min(len(n0), len(nb)); n0, nb = n0[:k], nb[:k]
    ch = np.loadtxt(nd + '/split/split_children.txt', dtype=int)[:k]
    act = nb[:, 7] > 0
    dang = np.full(k, np.nan); dang[act] = ang(n0[act, :3], nb[act, :3])
    dsh = np.hypot(nb[:, 3] - n0[:, 3], nb[:, 4] - n0[:, 4])
    corr, df, eo = nb[:, 8], 0.5 * (n0[:, 5] + n0[:, 6]), nb[:, 9]
    vb = views(n0[:, :3])
    cq = np.quantile(corr[act], [0.25, 0.5, 0.75])
    out.write('## %s: %d particles, %d refined; node corr quartiles %.4f %.4f %.4f\n' % (title, k, act.sum(), *cq))
    out.write('group      n   >30  >90  >150  med_ang  med_shift | corr p10 p50 p90 | defocus | eo odd | view TV vs rest | turned>30 share by corr quartile Q1..Q4\n')
    for g in np.unique(ch[ch > 0]):
        m = (ch == g) & act; o = (ch != g) & (ch > 0) & act
        a = dang[m]
        hv = np.bincount(vb[m], minlength=40) / m.sum(); ho = np.bincount(vb[o], minlength=40) / max(o.sum(), 1)
        q = np.digitize(corr[m], cq)
        byq = [np.mean(a[q == i] > 30) if (q == i).sum() else float('nan') for i in range(4)]
        out.write('%-6d %7d  %.3f %.3f %.3f  %6.2f  %6.2f | %.4f %.4f %.4f | %.3f | %.3f | %.3f | %s\n' % (
            g, m.sum(), np.mean(a > 30), np.mean(a > 90), np.mean(a > 150), np.median(a), np.median(dsh[m]),
            *np.quantile(corr[m], [0.1, 0.5, 0.9]), df[m].mean(), eo[m].mean(), 0.5 * np.abs(hv - ho).sum(),
            ' '.join('%.3f' % x for x in byq)))
    # the turned particles themselves: where do they sit in corr and in angle
    t = act & (dang > 30)
    out.write('turned > 30 in the node: %d (%.3f); their corr p50 %.4f vs untouched %.4f; angle histogram 30-60 %.3f 60-90 %.3f 90-150 %.3f 150-180 %.3f\n\n' % (
        t.sum(), t.mean(), np.median(corr[t]), np.median(corr[act & ~t]),
        *[np.mean((dang[t] >= lo) & (dang[t] < hi)) for lo, hi in ((30, 60), (60, 90), (90, 150), (150, 181))]))


T = sys.argv[1] + '/1_hetero_tree' if os.path.isdir(sys.argv[1] + '/1_hetero_tree') else sys.argv[1]
with open(sys.argv[2], 'w') as out:
    for nd, title in (('node_r', 'root r: groups 1..4 = r1 r2 r3 r4'), ('node_r4', 'node r4: groups 1..2 = r41 r42'), ('node_r1', 'node r1 (closed lineage): groups 1..3')):
        if os.path.exists(T + '/' + nd + '/split/split_children.txt'):
            profile(T + '/' + nd, out, title)
print('done')

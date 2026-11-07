#!/usr/bin/env python3
"""orc_prep.py OUT NODE DUMP TRUTH ALLTXT: inputs of the oracle refine=state_eval arms on a gate root (a project holding rows.txt's
rows; DUMP = print_project_field of its stage-B ptcl3D). Writes to OUT:
  lab_stageB.txt, lab_truth.txt, lab_coarse.txt   one state per project row (0 = inactive); truth = one label per truth class
                                                  (truth_map.txt), coarse = each class merged into its majority stage-B state
  prior_{stageB,truth,coarse}.txt                 SIMPLE_STATE_PRIOR rows 'iptcl lp_1..lp_K' (0 on the label, -1e4 elsewhere)
  poses_true.txt                                  e1 e2 e3 x y: the simulation's poses carried into the stage-B frame by the
                                                  global rotation and 3D shift fitted on the non-junk particles
  ptcl_info.txt                                   per row: combo junk jtype snrq bfac stageB_state ang_err shift_err (stage B vs truth)
  prep_report.txt                                 the frame fit and stage B's pose error."""
import os
import sys
import numpy as np

out, node, fdump, ftruth, fall = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
rep = open(os.path.join(out, 'prep_report.txt'), 'w')


def P(s=''):
    rep.write(s + '\n'); print(s, flush=True)


def read_kv(fn, keys):
    cols = {k: [] for k in keys}
    for ln in open(fn):
        if 'e1=' not in ln:
            continue
        d = dict(t.split('=', 1) for t in ln.split() if '=' in t)
        for k in keys:
            cols[k].append(float(d.get(k, 0.)))
    return {k: np.array(v) for k, v in cols.items()}


def rotmat(ang, axis):   # simple_ori_utils::rotmat (passive)
    c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
    r = np.zeros((len(ang), 3, 3))
    if axis == 2:
        r[:, 0, 0] = c; r[:, 0, 2] = -s; r[:, 1, 1] = 1.; r[:, 2, 0] = s; r[:, 2, 2] = c
    else:
        r[:, 0, 0] = c; r[:, 0, 1] = s; r[:, 1, 0] = -s; r[:, 1, 1] = c; r[:, 2, 2] = 1.
    return r


def euler2m(e):          # simple_ori_utils::euler2m: Rz(e3) Ry(e2) Rz(e1)
    return rotmat(e[:, 2], 3) @ rotmat(e[:, 1], 2) @ rotmat(e[:, 0], 3)


def m2euler(R):          # simple_ori_utils::m2euler
    c = np.clip(R[:, 2, 2], -1., 1.)
    e = np.zeros((len(R), 3))
    mid = (c < 0.9999999) & (c > -0.9999999); lo = c <= -0.9999999; hi = c >= 0.9999999
    e[mid, 0] = np.degrees(np.arctan2(R[mid, 2, 1], R[mid, 2, 0]))
    e[mid, 1] = np.degrees(np.arccos(c[mid]))
    e[mid, 2] = np.degrees(np.arctan2(R[mid, 1, 2], -R[mid, 0, 2]))
    e[lo, 0] = -np.degrees(np.arctan2(R[lo, 0, 1], R[lo, 1, 1])); e[lo, 1] = 180.
    e[hi, 0] = np.degrees(np.arctan2(R[hi, 0, 1], R[hi, 1, 1]))
    e[:, 0] %= 360.; e[:, 2] %= 360.
    return e


def geod(Ra, Rb):        # angle of Ra^T Rb, degrees
    E = np.einsum('nji,njk->nik', Ra, Rb)
    return np.degrees(np.arccos(np.clip((np.trace(E, axis1=1, axis2=2) - 1.) / 2., -1., 1.)))


def so3(M):
    U, _, Vt = np.linalg.svd(M)
    return U @ np.diag([1., 1., np.sign(np.linalg.det(U @ Vt))]) @ Vt


# the Euler conversions round-trip
rng = np.random.default_rng(0)
et = np.c_[rng.uniform(0, 360, 1000), np.degrees(np.arccos(rng.uniform(-1, 1, 1000))), rng.uniform(0, 360, 1000)]
assert np.abs(euler2m(m2euler(euler2m(et))) - euler2m(et)).max() < 1e-9

rows = np.loadtxt(os.path.join(node, 'rows.txt'), dtype=int, ndmin=1)
d = read_kv(fdump, ('e1', 'e2', 'e3', 'x', 'y', 'state', 'eo', 'indstk'))
n = len(d['e1'])
assert n == len(rows), (n, len(rows))
hdr = open(ftruth).readline().lstrip('#').split()
t = np.loadtxt(ftruth, ndmin=2)
assert (t[:, 0] == np.arange(1, len(t) + 1)).all(), 'truth rows are not 1..N'
T = t[rows - 1]
combo = T[:, hdr.index('combo')].astype(int); junk = T[:, hdr.index('junk')].astype(int) > 0
jtype = T[:, hdr.index('jtype')].astype(int); snrq = T[:, hdr.index('snrq')]; bfac = T[:, hdr.index('bfac')]
a = read_kv(fall, ('e1', 'e2', 'e3', 'x', 'y'))
assert len(a['e1']) >= rows.max()
et = np.c_[a['e1'], a['e2'], a['e3']][rows - 1]; xt = np.c_[a['x'], a['y']][rows - 1]
eb = np.c_[d['e1'], d['e2'], d['e3']]; xb = np.c_[d['x'], d['y']]
sB = d['state'].astype(int); act = sB > 0
P('node %s: %d rows, %d active, stage-B states %s; true shifts |x|,|y| max %.3f px' % (node, n, act.sum(), np.unique(sB[act]).tolist(), np.abs(xt).max()))
if not (np.unique(d['indstk']) .size == n):
    P('WARNING: indstk is not one-to-one over the rows')

Rb, Rt = euler2m(eb), euler2m(et)
raw = geod(Rb, Rt)
good = act & ~junk & (raw < 15.)
fits = {}
for side in ('left', 'right'):
    D = np.einsum('nij,nkj->nik', Rb, Rt) if side == 'left' else np.einsum('nji,njk->nik', Rt, Rb)
    G = so3(D[good].mean(0))
    for _ in range(3):   # refit on the particles near the fitted offset
        dev = geod(np.broadcast_to(G, D.shape), D)
        sel = good & (dev < 10.)
        G = so3(D[sel].mean(0))
    dev = geod(np.broadcast_to(G, D.shape), D)
    fits[side] = (np.median(dev[good]), G, geod(G[None], np.eye(3)[None])[0])
    P('global rotation (%s): |G| %.3f deg, median residual angle %.3f deg (raw %.3f)' % (side, fits[side][2], fits[side][0], np.median(raw[good])))
side = min(fits, key=lambda s: fits[s][0])
G = fits[side][1]
RtB = G @ Rt if side == 'left' else Rt @ G
# 3D translation of the frame: the residual shifts as the projection of one vector, through R or R^T
r = xb - xt
best = None
for model in ('R', 'RT'):
    A = RtB[:, :2, :] if model == 'R' else np.transpose(RtB, (0, 2, 1))[:, :2, :]
    gs = good & (np.hypot(r[:, 0], r[:, 1]) < 6.)
    tv = np.linalg.solve(np.einsum('nij,nik->jk', A[gs], A[gs]), np.einsum('nij,ni->j', A[gs], r[gs]))
    res = r - np.einsum('nij,j->ni', A, tv)
    rms0, rms1 = np.sqrt((r[gs] ** 2).sum(1).mean()), np.sqrt((res[gs] ** 2).sum(1).mean())
    P('frame shift (%s): t = [%.3f %.3f %.3f] px, rms shift error %.3f -> %.3f px' % (model, *tv, rms0, rms1))
    if best is None or rms1 < best[0]:
        best = (rms1, model, tv, A)
_, model, tv, A = best
eB = m2euler(RtB)
xB = xt + np.einsum('nij,j->ni', A, tv)
P('true poses in the stage-B frame: rotation %s (%.3f deg), translation model %s |t| %.3f px' % (side, fits[side][2], model, np.linalg.norm(tv)))
np.savetxt(os.path.join(out, 'poses_true.txt'), np.c_[eB, xB], fmt='%.6f')
ang = geod(Rb, euler2m(eB)); sh = np.hypot(*(xb - xB).T)


def stats(m):
    return 'n %6d  angle median %5.2f mean %5.2f p90 %5.2f  >5 %.3f >10 %.3f >30 %.3f | shift median %.2f p90 %.2f px' % (
        m.sum(), np.median(ang[m]), ang[m].mean(), np.percentile(ang[m], 90), (ang[m] > 5).mean(), (ang[m] > 10).mean(),
        (ang[m] > 30).mean(), np.median(sh[m]), np.percentile(sh[m], 90))


P('stage-B pose error against the truth (in the stage-B frame):')
P('  all        ' + stats(act))
P('  non-junk   ' + stats(act & ~junk))
P('  junk       ' + stats(act & junk))
for c in np.unique(combo):
    P('  class %3d  ' % c + stats(act & (combo == c)))

# labels; ORC_MINN / ORC_MINFRAC (a child node): classes below max(MINN, MINFRAC x active) join junk (0) as one 'other' reference
minn = max(float(os.environ.get('ORC_MINN', 0)), float(os.environ.get('ORC_MINFRAC', 0)) * act.sum())
cnt = {c: (combo[act] == c).sum() for c in np.unique(combo[act])}
rare = sorted(int(c) for c in cnt if c != 0 and cnt[c] < minn)
combo_ref = np.where(np.isin(combo, rare), 0, combo)
if rare:
    P('classes below %d particles, referenced with junk as one other class: %s' % (minn, {c: int(cnt[c]) for c in rare}))
classes = np.unique(combo_ref[act])
lab_t = np.zeros(n, int); lab_t[act] = np.searchsorted(classes, combo_ref[act]) + 1
maj = {c: int(np.bincount(sB[act & (combo_ref == c)]).argmax()) for c in classes}
used = sorted(set(maj.values()))
lab_c = np.zeros(n, int); lab_c[act] = [used.index(maj[c]) + 1 for c in combo_ref[act]]
with open(os.path.join(out, 'truth_map.txt'), 'w') as f:
    for k, c in enumerate(classes):
        f.write('%d %d %d\n' % (k + 1, c, (combo_ref[act] == c).sum()))
with open(os.path.join(out, 'coarse_map.txt'), 'w') as f:
    for c in classes:
        f.write('class %d -> stage-B state %d -> coarse %d\n' % (c, maj[c], used.index(maj[c]) + 1))
P('truth classes %s -> labels 1..%d; coarse: %s (%d groups)' % (classes.tolist(), len(classes), {int(c): maj[c] for c in classes}, len(used)))
for name, lab in (('stageB', sB), ('truth', lab_t), ('coarse', lab_c)):
    np.savetxt(os.path.join(out, 'lab_%s.txt' % name), lab, fmt='%d')
    K = lab.max()
    with open(os.path.join(out, 'prior_%s.txt' % name), 'w') as f:
        for i in np.nonzero(lab > 0)[0]:
            f.write('%d %s\n' % (i + 1, ' '.join('0' if k == lab[i] else '-1e4' for k in range(1, K + 1))))
    P('labels %-7s K %d sizes %s' % (name, K, np.bincount(lab[act])[1:].tolist()))
np.savetxt(os.path.join(out, 'ptcl_info.txt'), np.c_[combo, junk, jtype, snrq, bfac, sB, ang, sh],
           fmt=['%d', '%d', '%d', '%.4f', '%.2f', '%d', '%.3f', '%.3f'],
           header='combo junk jtype snrq bfac stageB_state ang_err_deg shift_err_px')
rep.close()

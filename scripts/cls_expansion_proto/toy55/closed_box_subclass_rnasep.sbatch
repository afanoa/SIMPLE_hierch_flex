#!/bin/bash
#SBATCH -J rnasep_closed_box_subcls
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Two diagnostics for "the helix is never long enough", python only. (A) Where the closed arm lives relative to our box:
# radial mass profile of the deposited-18 arm in our frame (box 234 A, inscribed radius 117 A), the tip metric of the
# deposited map itself after SIMPLE-style soft masks (200/220/232/240 A) and the inscribed sphere; and the lever-arm
# blur: the deposited map averaged over random pose errors (rotation vectors with median angle 2/4/6 deg, shifts median
# 1.4 px) -> tip occupancy and contact visibility expected from pose errors alone. (B) Sub-class pooling of the
# per-particle focal score on the family-1 board: how much of its variance sits between the 2D expansion sub-classes
# (vs a permutation null), the sub-classes the 3D score calls most docked (montage of their averages for the eye), and
# a re-ranking with the shrunk sub-class mean as prior (states_subcls.txt for a later reconstruction).
# env: FD (board dir), N1 (30000)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; FD=${FD:-$KF/closed_focal}; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs.mrc; [ -s $CAVGS ] || CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs_even.mrc
D=$FD/box_subcls; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') bx] $*"; }; export OMP_NUM_THREADS=8; N1=${N1:-30000}
log "(A) the closed arm vs the box, and the lever-arm blur of the deposited map"
$PY - <<PYEOF
import numpy as np, glob, sys, time
sys.path.insert(0, '$W'); from model_class import read_mrc
from scipy.ndimage import gaussian_filter, map_coordinates
W = '$W'; SMPD = 1.464; n = 160; rng = np.random.default_rng(0)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0])
g = lambda v: gaussian_filter(v, 2.5); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std())
zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2); r = rr * SMPD; tip = arm & (r > 95); base = arm & (r < 80)
print('arm %d voxels, tip (r > 95 A) %d, core %d; arm radius: median %.0f A, 90th pct %.0f, max %.0f; dense closed map max radius %.0f A' % (arm.sum(), tip.sum(), core.sum(), np.median(r[arm]), np.percentile(r[arm], 90), r[arm].max(), r[dens18].max()))
edges = [0, 60, 70, 80, 90, 95, 100, 105, 110, 114, 117, 130]; tot = L18[arm].sum()
print('arm mass by radius: ' + ', '.join('%d-%d A %.1f %%' % (lo, hi, 100 * L18[arm & (r >= lo) & (r < hi)].sum() / tot) for lo, hi in zip(edges[:-1], edges[1:])))
rt = L18[tip].mean() / L18[core].mean(); ra = L18[arm].mean() / L18[core].mean()
def tipocc(v): Lv = g(v); c = Lv[core].mean(); return (Lv[tip].mean() / c) / rt, (Lv[arm].mean() / c) / ra
lp = lambda v, res: np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v))) * (rr <= n * SMPD / res)))))
def vis(v): L = lp(v, 8.); L = L / max(np.median(L[core]), 1e-30); return np.percentile(L[tip], 99), np.percentile(L[tip], 90), np.mean(L[tip] > 0.5)
def softmask(v, diam, edge_px=6.): rad = diam / 2. / SMPD; w = np.clip((rad - rr) / edge_px, 0, 1); return v * (0.5 * (1 - np.cos(np.pi * w)))
print('deposited 18 itself: tip %.2f arm %.2f | visibility p99 %.2f p90 %.2f >0.5 %.2f' % (*tipocc(c18), *vis(c18)))
for d in (200, 220, 232, 240): t, a_ = tipocc(softmask(c18, d)); print('  after a SIMPLE-style soft mask of %d A (6-px cosine edge): tip %.2f arm %.2f | visibility p99 %.2f p90 %.2f >0.5 %.2f' % (d, t, a_, *vis(softmask(c18, d))))
t, a_ = tipocc(c18 * (r <= 117)); print('  after a hard inscribed sphere (117 A, what an in-plane rotation of the box keeps): tip %.2f arm %.2f | visibility p99 %.2f p90 %.2f >0.5 %.2f' % (t, a_, *vis(c18 * (r <= 117))))
coords = np.stack([zz - n / 2, yy - n / 2, xx - n / 2], 0).reshape(3, -1).astype(np.float32)
def rotvec2m(v):
    a = np.linalg.norm(v)
    if a < 1e-12: return np.eye(3)
    k = v / a; K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]]); return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K
def blur(v, med_deg, med_shift_px, ns=120):
    sig = np.radians(med_deg) / 1.538; ssig = med_shift_px / 1.177; acc = np.zeros_like(v); angs = []
    for _ in range(ns):
        rv = rng.normal(size=3) * sig; R = rotvec2m(rv); s = rng.normal(size=3) * ssig * np.sqrt(2. / 3.); angs.append(np.degrees(np.linalg.norm(rv)))
        c = (R @ coords) + s[:, None] + n / 2; acc += map_coordinates(v, c.reshape(3, n, n, n), order=1, mode='constant', cval=0.)
    return acc / ns, np.median(angs), np.percentile(angs, 90)
print('lever-arm blur of the deposited closed map (isotropic pose errors, 120 samples; toy pipeline: angle median 3.98 deg, p90 7.35, shift median 1.44 px):')
for md in (2., 4., 6.):
    t0 = time.time(); bv, m_, p90 = blur(c18, md, 1.44); t, a_ = tipocc(bv); print('  angle median %.1f deg (p90 %.1f), shift median 1.4 px: tip %.2f arm %.2f | visibility p99 %.2f p90 %.2f >0.5 %.2f   (%.0f s)' % (m_, p90, t, a_, *vis(bv), time.time() - t0))
bv, m_, p90 = blur(c18, 4., 0.); t, a_ = tipocc(bv); print('  angle median %.1f deg, no shift error: tip %.2f arm %.2f | visibility p99 %.2f p90 %.2f >0.5 %.2f' % (m_, t, a_, *vis(bv)))
bv, m_, p90 = blur(c18, 0.01, 1.44); t, a_ = tipocc(bv); print('  shift error only (median 1.4 px): tip %.2f arm %.2f | visibility p99 %.2f p90 %.2f >0.5 %.2f' % (t, a_, *vis(bv)))
PYEOF
log "(B) sub-class pooling of the per-particle focal score on the board $FD"
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M = np.load('$META'); flags = M['flags'].astype(int); lab = M['lab692'].astype(int); N = len(flags); ref = np.where(flags == 1)[0]
F = np.loadtxt('$FD/board/focal_scores.txt', comments='#'); fs = F[:, 1]; R = np.load('$FD/board/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(R['i1']); s = r[:, :ngrid] @ fs[:ngrid]
pos = -np.ones(N, int); pos[ref] = np.arange(len(ref)) if len(lab) == len(ref) else np.arange(len(ref)); sl = lab[pos[pinds - 1]] if len(lab) == len(ref) else lab[pinds - 1]
ok = (pos[pinds - 1] >= 0) if len(lab) == len(ref) else np.ones(len(s), bool); s = s[ok]; sl = sl[ok]; pinds = pinds[ok]; print('board particles with a sub-class label: %d (labels array %d, flags==1 rows %d)' % (ok.sum(), len(lab), len(ref)))
u, inv = np.unique(sl, return_inverse=True); cnt = np.bincount(inv); K = len(u); mu = s.mean(); mean = np.bincount(inv, s) / cnt
tot = s.var(); between = np.sum(cnt * (mean - mu) ** 2) / len(s); within = tot - between
nulls = []
for _ in range(20):
    sp = rng = np.random.default_rng(_).permutation(s); mp = np.bincount(inv, sp) / cnt; nulls.append(np.sum(cnt * (mp - mu) ** 2) / len(s))
print('focal score variance: total %.4f, between %d sub-classes %.4f (%.1f %%), permutation null %.4f (%.1f %%) -> excess %.1f %% of the variance is sub-class structure' % (tot, K, between, 100 * between / tot, np.mean(nulls), 100 * np.mean(nulls) / tot, 100 * (between - np.mean(nulls)) / tot))
tau2 = max(between - np.mean(nulls), 1e-6); shrunk = (cnt * mean / within + mu / tau2) / (cnt / within + 1. / tau2)
corr_self = np.corrcoef(s, shrunk[inv])[0, 1]; print('correlation of a particle score with its (shrunk) sub-class mean: %.3f; sub-class mean sd %.3f vs particle score sd %.3f' % (corr_self, shrunk.std(), s.std()))
o = np.argsort(-shrunk); print('most docked sub-classes by shrunk mean focal score (sub-class, members on the board, mean):'); print('  ' + ' '.join('%d(n%d,%+.2f)' % (u[k], cnt[k], shrunk[k]) for k in o[:16])); print('least docked:'); print('  ' + ' '.join('%d(n%d,%+.2f)' % (u[k], cnt[k], shrunk[k]) for k in o[-8:]))
N1 = $N1; comb = s + shrunk[inv]; top_s = set(pinds[np.argsort(-s)[:N1]]); top_c = set(pinds[np.argsort(-comb)[:N1]]); print('top %d by own score vs by own + sub-class prior: overlap %d (%.1f %%)' % (N1, len(top_s & top_c), 100. * len(top_s & top_c) / N1))
st = np.zeros(N, int); st[np.array(sorted(top_c)) - 1] = 1; np.savetxt('$D/states_subcls.txt', st, fmt='%d'); np.savetxt('$D/subclass_focal_means.txt', np.c_[u, cnt, mean, shrunk], fmt='%d %d %.4f %.4f', header='subclass n_on_board mean_focal shrunk_mean')
# montage of the 16 most and 8 least docked sub-class averages, for the eye
try:
    C = read_mrc('$CAVGS'); nav = C.shape[0]; sel = [u[k] for k in o[:16]] + [u[k] for k in o[-8:]]; fig, ax = plt.subplots(3, 8, figsize=(16, 6.4))
    for i, a in enumerate(ax.ravel()):
        k = sel[i]; img = C[k - 1] if k - 1 < nav else np.zeros((C.shape[1], C.shape[2])); a.imshow(img, cmap='gray'); a.set_title('%s sub %d n%d %+.2f' % ('docked' if i < 16 else 'undocked', k, cnt[np.where(u == k)[0][0]], shrunk[np.where(u == k)[0][0]]), fontsize=7); a.axis('off')
    plt.suptitle('2D expansion sub-class averages the 3D focal score calls most docked (rows 1-2) and least docked (row 3)'); plt.tight_layout(); plt.savefig('$D/subclass_focal_montage.png', dpi=110); print('montage -> $D/subclass_focal_montage.png (sub-class k = average k, 1-based)')
except Exception as e: print('montage skipped: %s' % e)
PYEOF
log DONE

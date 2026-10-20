#!/usr/bin/env python3
"""Scale-normalise the docked deposited maps (apo + holo) for the apo/holo toy: EMDB maps carry arbitrary scales, and a
holo/apo scale difference would hand a classifier a trivial brightness axis. Every map is scaled by least squares onto
the apo consensus over the RIGID CORE (consensus density above 1 sd, outside the moving arm and outside the holo excess),
written to OUT_DIR under its original name. Then the protein site: mean (holo - apo consensus) after scaling, positive,
outside the arm, above 3 sd, largest connected component, dilated and softened -> rnpa_mask.mrc (+ the mean excess map),
with its volume and radius reported. usage: normalize_ah_maps.py APO_DOCKED HOLO_DOCKED CONSENSUS OUT_DIR MASK_OUT_DIR"""
import sys, os, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation, label
apo_dir, holo_dir, cons_fn, out, mout = sys.argv[1:6]; os.makedirs(out, exist_ok=True); os.makedirs(mout, exist_ok=True); SMPD = 1.464
cons = read_mrc(cons_fn); g = lambda v: gaussian_filter(v, 2.0); Lc = g(cons); dens = Lc > Lc.mean() + 1.0 * Lc.std()
apo = {os.path.basename(f): read_mrc(f) for f in sorted(glob.glob(apo_dir + '/conf*_m160.mrc'))}; holo = {os.path.basename(f): read_mrc(f) for f in sorted(glob.glob(holo_dir + '/conf*_m160.mrc'))}
# the moving arm: across-apo variance; the holo excess (raw): both excluded from the scaling core
La = {k: g(v) for k, v in apo.items()}; mean_a = np.mean(list(La.values()), 0); var_a = np.mean([(v - mean_a) ** 2 for v in La.values()], 0)
arm = var_a > np.percentile(var_a[dens], 80); arm = binary_dilation(arm & dens, iterations=3)
def scale_to(v, core):
    x = g(v)[core]; y = Lc[core]; return float(np.dot(x, y) / max(np.dot(x, x), 1e-30))
core0 = dens & ~arm
raw_ex = np.mean([np.maximum(g(v) * scale_to(v, core0) - Lc, 0) for v in holo.values()], 0); exc = raw_ex > raw_ex.mean() + 3 * raw_ex.std(); exc = binary_dilation(exc, iterations=3)
core = dens & ~arm & ~exc; print('scaling core: %d voxels (density %d, arm %d, raw holo excess %d excluded)' % (core.sum(), dens.sum(), (arm & dens).sum(), (exc & dens).sum()))
sc = {}
for name, v in list(apo.items()) + list(holo.items()):
    s = scale_to(v, core); sc[name] = s; write_mrc(out + '/' + name, (v * s).astype(np.float32), SMPD)
print('scale factors onto the consensus core: apo %s | holo %s' % (' '.join('%s:%.2f' % (k.split('_')[0], sc[k]) for k in apo), ' '.join('%s:%.2f' % (k.split('_')[0], sc[k]) for k in holo)))
E = np.mean([np.maximum(g(v) * sc[k] - Lc, 0) for k, v in holo.items()], 0); E[arm] = 0.
thr = E.mean() + 3 * E.std(); m = E > thr; lab, nl = label(m); sizes = np.bincount(lab.ravel())[1:]
if nl == 0: raise SystemExit('no holo excess above 3 sd')
order = np.argsort(-sizes); n = cons.shape[0]; zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD
for j in order[:4]:
    c = lab == j + 1; print('  excess component %d: %d voxels = %.0f A^3, mass share %.2f, centroid radius %.0f A, max radius %.0f A' % (j + 1, c.sum(), c.sum() * SMPD ** 3, E[c].sum() / E[m].sum(), np.mean(rr[c]), rr[c].max()))
site = binary_dilation(lab == order[0] + 1, iterations=2); soft = gaussian_filter(site.astype(float), 2.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
write_mrc(mout + '/rnpa_mask.mrc', soft.astype(np.float32), SMPD); write_mrc(mout + '/holo_mean_excess.mrc', E.astype(np.float32), SMPD); write_mrc(mout + '/arm_mask.mrc', arm.astype(np.float32), SMPD)
print('rnpa_mask.mrc = largest excess component dilated: %d voxels above 0.5 = %.0f A^3 (rnpA ~ 16,000 A^3 at 1.35 g/cm3), centroid radius %.0f A' % ((soft > 0.5).sum(), (soft > 0.5).sum() * SMPD ** 3, np.mean(rr[soft > 0.5])))
# per holo map: how much of the site each map fills (occupancy readout of the deposited maps themselves)
for k, v in holo.items(): d = g(v) * sc[k] - Lc; print('  %s: site density relative to the mean holo excess %.2f' % (k.split('_')[0], float(np.sum(soft * np.maximum(d, 0)) / max(np.sum(soft * E), 1e-30))))

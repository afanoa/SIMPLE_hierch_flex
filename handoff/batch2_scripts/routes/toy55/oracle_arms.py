#!/usr/bin/env python3
"""Per-particle ORACLE ceiling of every arm angle of the IgM-like set from arm_llr.py's scores of each image against ALL true conformer maps
(arm_llr.npz: SSE per image x conformer, whitened, per-image LS scale): posterior over the conformers p(k | x) ~ exp(-SSE / tau) (uniform
prior; tau scanned, the best per angle reported = an optimistic calibration), posterior-mean angle per arm, and its correlation with the
true angle (|corr| and R^2), whole image and masked; plus the MAP conformer's angle. Also the binary ceiling for the sign of each angle
(AUC of P(angle > 0 | x)). usage: oracle_arms.py ARM_LLR.npz [ARM_LLR2.npz ...]"""
import sys, numpy as np
for fn in sys.argv[1:]:
    z = np.load(fn, allow_pickle=True); TH, kt = z['TH'], z['ktrue']; narm = TH.shape[1] // 2; truth = TH[kt]
    print('== %s: %d images, %d conformers, %d arms' % (fn, len(kt), TH.shape[0], narm))
    for key in ('SSE', 'SSEm'):
        S = z[key].astype(np.float64); best = {}
        for tau in (0.25, 0.5, 1., 2., 4., 8.):
            L = -S / tau; L -= L.max(1, keepdims=True); P = np.exp(L); P /= P.sum(1, keepdims=True); E = P @ TH
            for j in range(TH.shape[1]):
                r = np.corrcoef(E[:, j], truth[:, j])[0, 1]
                if j not in best or r > best[j][0]: best[j] = (r, tau)
        mp = TH[np.argmin(S, 1)]; line = []
        for k in range(narm):
            for t, nm in ((0, 'th'), (1, 'ph')):
                j = 2 * k + t; rm = np.corrcoef(mp[:, j], truth[:, j])[0, 1]; line.append('a%d_%s %.2f (tau %g, MAP %.2f)' % (k + 1, nm, best[j][0], best[j][1], rm))
        print('  %s posterior-mean |corr| with the true angle: %s' % ('whole image' if key == 'SSE' else 'masked to the arm region', '; '.join(line)))

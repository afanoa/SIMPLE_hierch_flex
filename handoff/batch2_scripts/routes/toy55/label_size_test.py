#!/usr/bin/env python3
"""usage: label_size_test.py TRUTH NAME:post_master:curve [..]: label-set size for a site pole from its hold-out curve: top-N* (argmax of the
hold-out correlation = the pole MAP's SNR optimum) vs top-EPM (the saturation of the equivalent pure members = the estimated member count);
truth precision / recall of each label set, and the truth-optimal F1 for reference."""
import sys, numpy as np
TJ = sys.argv[1]; hd = [l.split()[1:] for l in open(TJ) if l.startswith('#')][0]; T = np.loadtxt(TJ); c = lambda n: T[:, hd.index(n)]
y = (c('C_closed') > 0) & (c('junk') == 0); M = y.sum()
for a in sys.argv[2:]:
    nm, pf, cf = a.split(':'); P = np.loadtxt(pf); C = np.loadtxt(cf); j = int(np.argmax(C[:, 2])); nstar = int(C[j, 0]); epm = int(round(C[:, 6].max()))
    o = np.argsort(-np.where(P >= 0, P, -np.inf)); cy = np.cumsum(y[o]); f1 = 2 * cy / (np.arange(1, len(o) + 1) + M); nb = int(np.argmax(f1)) + 1
    def pr(N): return cy[N - 1] / N, cy[N - 1] / M, 2 * cy[N - 1] / (N + M)
    print('%-12s N* %6d: precision %.3f recall %.3f F1 %.3f | EPM %6d: precision %.3f recall %.3f F1 %.3f | truth-best N %6d: F1 %.3f (true members %d)'
          % ((nm, nstar) + pr(nstar) + (epm,) + pr(epm) + (nb, f1[nb - 1], M)))

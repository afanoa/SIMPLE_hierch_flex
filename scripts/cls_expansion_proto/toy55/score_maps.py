#!/usr/bin/env python3
"""Score a set of full-box real-data maps against the 10 docked deposited conformers the way deposited_match.py scores
cells: the deviation of every map from the MEAN of the set (the common real-vs-toy offset cancels) against every
conformer's deviation from the consensus, correlation in a Fourier band; per map: own-conformer corr (if the maps are
labelled by conformer), best conformer, rank of the own one. Bands 60-20 A and 30-12 A.
usage: score_maps.py MAPS_DIR [--labels 1,10,18,...] [--pattern 'recvol_state%02d.mrc'] [--n 10]"""
import sys, glob, re, argparse, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, Shells, fftc
ap = argparse.ArgumentParser(); ap.add_argument('maps'); ap.add_argument('--labels', default=''); ap.add_argument('--pattern', default='recvol_state%02d.mrc'); ap.add_argument('--n', type=int, default=10); a = ap.parse_args()
W = '/mnt/beegfs/elmlund/afan/routes/toy55'; SMPD = 1.464; ids = [1, 10, 18, 27, 28, 37, 46, 48, 49, 51]; S = Shells(160, SMPD); msk = S.soft_mask(100 / SMPD)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in ids}; Cm = np.mean([C[i] for i in ids], 0)
V = [read_mrc('%s/%s' % (a.maps, a.pattern % (k + 1))) for k in range(a.n)]; Vm = np.mean(V, 0); labels = [int(x) for x in a.labels.split(',')] if a.labels else [None] * a.n
def bc(x, y, lo, hi):
    m = S.band(lo, hi); return float(np.real(np.sum(np.conj(x[m]) * y[m])) / max(np.sqrt(np.sum(np.abs(x[m]) ** 2) * np.sum(np.abs(y[m]) ** 2)), 1e-30))
for ref_name, ref in (('conformer - consensus', cons), ('conformer - mean of conformers', Cm)):
    T = {i: fftc((C[i] - ref) * msk) for i in ids}
    for lo, hi in ((60, 20), (30, 12)):
        print('== maps - mean of maps vs %s, %g-%g A' % (ref_name, lo, hi)); hits = 0; own_all = []
        for k in range(a.n):
            D = fftc((V[k] - Vm) * msk); cc = {i: bc(D, T[i], lo, hi) for i in ids}; order = sorted(ids, key=lambda i: -cc[i]); best = order[0]; lab = labels[k]
            own = cc[lab] if lab else float('nan'); rank = order.index(lab) + 1 if lab else 0; hits += (lab == best); own_all.append(own)
            print('  map %2d%s: best %2d %+.2f%s | %s' % (k + 1, (' (conformer %2d)' % lab) if lab else '', best, cc[best], (' | own %+.2f rank %d' % (own, rank)) if lab else '', ' '.join('%d:%+.2f' % (i, cc[i]) for i in ids)))
        if labels[0]: print('  own conformer is the best match for %d of %d maps; mean own corr %+.2f' % (hits, a.n, np.nanmean(own_all)))

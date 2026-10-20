#!/usr/bin/env python3
"""Candidates for validating the f692e states on the independent halves: every half row gets its f692e leaf label
(consolidated_12_fullw_region, 63 states), except the zoomed closed sub-states, which get their own labels (pure closed
zoom state 1 -> CLOSED_PURE; zoom states 143, 77, 30 -> CLOSED_NEAR). Leaves with fewer than MINC particles in either half
become 0 and are left out of the validation (the bench loads labelled rows only). usage: f692e_candidates.py REGION_STATES ZOOM_STATES SIM_PREFIX OUT_PREFIX [MINC]"""
import sys, numpy as np
regf, zoomf, simp, outp = sys.argv[1:5]; minc = int(sys.argv[5]) if len(sys.argv) > 5 else 1500
R = np.loadtxt(regf).astype(int); Z = np.loadtxt(zoomf).astype(int); PURE, NEAR = 64, 65
lab = R.copy(); lab[Z == 1] = PURE; lab[np.isin(Z, [143, 77, 30])] = NEAR
print('meta rows %d: f692e leaves on %d rows; closed pure %d, closed near %d (taken out of their leaves: %d / %d rows had a leaf label)' % (len(R), (R > 0).sum(), (lab == PURE).sum(), (lab == NEAR).sum(), ((Z == 1) & (R > 0)).sum(), (np.isin(Z, [143, 77, 30]) & (R > 0)).sum()))
H = {h: lab[np.loadtxt('%s_%s/meta_rows.txt' % (simp, h)).astype(int) - 1] for h in 'AB'}
cnt = {h: np.bincount(H[h], minlength=66) for h in 'AB'}; small = [k for k in range(1, 64) if min(cnt['A'][k], cnt['B'][k]) < minc]
for h in 'AB':
    L = H[h].copy(); L[np.isin(L, small)] = 0; np.savetxt('%s_%s.txt' % (outp, h), L, fmt='%d')
    print('half %s: %d rows, %d candidates in %d groups (%d leaves below %d left out), closed pure %d, closed near %d' % (h, len(L), (L > 0).sum(), len(set(L[L > 0])), len(small), minc, (L == PURE).sum(), (L == NEAR).sum()))
print('leaves below the floor: %s' % ' '.join('%d(%d/%d)' % (k, cnt['A'][k], cnt['B'][k]) for k in small))

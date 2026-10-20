#!/usr/bin/env python3
"""Borrow the core from the parent: a state's map outside its own region is the parent's map (every particle of the parent
node, i.e. the full-SNR core), inside the region it is the state's own reconstruction. Composed per half:
M_h = P_h + m (S_h - P_h), with S scale-matched to P on the support and m the state's region mask (soft, 15 A).
Readouts per state: global FSC of the composed map, FSC inside the region (state-only vs composed), both at 0.143.
usage: compose_state_maps.py RUN CONSOLIDATED_DIR PARENT_PREFIX OUT [--smpd 1.464]   (PARENT_PREFIX + _even.mrc / _odd.mrc / .mrc)"""
import sys, os, re, glob, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc, write_mrc
ap = argparse.ArgumentParser(); ap.add_argument('run'); ap.add_argument('cons'); ap.add_argument('parent'); ap.add_argument('out'); ap.add_argument('--smpd', type=float, default=1.464); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
PAR = {}
def parent_of(fam):   # per-family parent when the prefix carries {F} (each family's states live in that family's refined frame)
    key = a.parent.replace('{F}', fam)
    if key not in PAR: PAR[key] = {h: read_mrc(key + s + '.mrc').astype(np.float32) for h, s in (('e', '_even'), ('o', '_odd'), ('f', ''))}
    return PAR[key]
P = parent_of('1'); n = P['f'].shape[0]
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2; sh = np.minimum(np.round(np.sqrt(k2) * n).astype(int), n // 2).ravel()
lp = lambda v, res: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (a.smpd / res) ** 2))).astype(np.float32)
Pl = lp(P['f'], 10.); supp = Pl > 0.2 * Pl.max(); x = np.arange(n) - n // 2; rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); sph = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1).astype(np.float32)
def fsc_res(e, o, m):
    E, O = np.fft.fftn(e * m), np.fft.fftn(o * m); num = np.bincount(sh, np.real(np.conj(E) * O).ravel(), n // 2 + 1); d1 = np.bincount(sh, (np.abs(E) ** 2).ravel(), n // 2 + 1); d2 = np.bincount(sh, (np.abs(O) ** 2).ravel(), n // 2 + 1)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * a.smpd / k
def scale_to_parent(S):
    A = np.stack([Pl[supp], np.ones(supp.sum())], 1); ab = np.linalg.lstsq(A, lp(S, 10.)[supp], rcond=None)[0]; return (S - ab[1]) / max(ab[0], 1e-12)
_mc = {}
def region_mask(leaf):
    m = re.match(r'f(\d)_(motion_)?_?(\d+)', leaf); fam, mot, reg = m.group(1), m.group(2), m.group(3); fn = '%s/fam%s/%s/region_%s.mrc' % (a.run, fam, 'regions_motion' if mot else 'regions', reg)
    if fn not in _mc: _mc[fn] = np.clip(lp(read_mrc(fn).astype(np.float32), 15.), 0, 1) if os.path.isfile(fn) else None
    return _mc[fn], fn
groups = {int(l.split()[0]): l.split()[2] for l in open(a.cons + '/groups.txt') if not l.startswith('#')}
print('state  particles  leaf  region_mask_fraction  FSC_global(state->composed)  FSC_in_region(state->composed)  FSC_core_outside(state->composed)')
for s in sorted(groups):
    leaf = groups[s]; mk, fn = region_mask(leaf)
    if mk is None: print('%d %s: region mask missing (%s)' % (s, leaf, fn)); continue
    P = parent_of(leaf[1]); Pl = lp(P['f'], 10.); supp = Pl > 0.2 * Pl.max()
    S = {h: read_mrc('%s/fam0/regionall/level1/maps/recvol_state%02d%s.mrc' % (a.cons, s, suf)).astype(np.float32) for h, suf in (('e', '_even'), ('o', '_odd'), ('f', ''))}
    if not all(os.path.isfile('%s/fam0/regionall/level1/maps/recvol_state%02d%s.mrc' % (a.cons, s, suf)) for suf in ('_even', '_odd')): continue
    S = {h: scale_to_parent(v) for h, v in S.items()}; M = {h: P[h] + mk * (S[h] - P[h]) for h in S}
    write_mrc('%s/state_%02d.mrc' % (a.out, s), M['f'].astype(np.float32), a.smpd); write_mrc('%s/state_%02d_even.mrc' % (a.out, s), M['e'].astype(np.float32), a.smpd); write_mrc('%s/state_%02d_odd.mrc' % (a.out, s), M['o'].astype(np.float32), a.smpd)
    npart = [int(l.split()[1]) for l in open(a.cons + '/groups.txt') if not l.startswith('#') and int(l.split()[0]) == s][0]
    print('%2d  %7d  %-18s %.2f   %.1f -> %.1f   %.1f -> %.1f   %.1f -> %.1f' % (s, npart, leaf, float(mk.mean()), fsc_res(S['e'], S['o'], sph), fsc_res(M['e'], M['o'], sph), fsc_res(S['e'], S['o'], mk ** 2 * sph), fsc_res(M['e'], M['o'], mk ** 2 * sph), fsc_res(S['e'], S['o'], (1 - mk) ** 2 * sph), fsc_res(M['e'], M['o'], (1 - mk) ** 2 * sph)), flush=True)
for key, Pk in PAR.items(): print('parent %s: global FSC %.1f A' % (key, fsc_res(Pk['e'], Pk['o'], sph)))

#!/bin/bash
#SBATCH -J klem_diag_mask
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# why does the region mask inflate the half-map FSC of a leaf? mask value range, spectral content, own FSC under mask variants
/mnt/nasapps/production/python/3.12/bin/python3 - <<'PYEOF'
import numpy as np, sys
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc
A = '/mnt/beegfs/elmlund/afan/rnasep_klem/f692c'; smpd = 1.464
cons = read_mrc('/mnt/beegfs/elmlund/afan/routes/toy55/maps160/toy_consensus.mrc').astype(np.float32); n = cons.shape[0]
f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2; sh = np.minimum(np.round(np.sqrt(k2) * n).astype(int), n // 2).ravel()
G = lambda res: np.exp(-0.5 * k2 * (res / 1.2) ** 2 / smpd ** 2 * smpd ** 2)   # as in merge_data_derived (res in A -> px via smpd... printed below)
lp = lambda v, res: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-0.5 * k2 * (res / smpd / 1.2) ** 2)))
def fsc(e, o):
    E, O = np.fft.fftn(e), np.fft.fftn(o); num = np.bincount(sh, np.real(np.conj(E) * O).ravel(), n // 2 + 1); d1 = np.bincount(sh, (np.abs(E) ** 2).ravel(), n // 2 + 1); d2 = np.bincount(sh, (np.abs(O) ** 2).ravel(), n // 2 + 1)
    return num / np.maximum(np.sqrt(d1 * d2), 1e-30)
def spec_frac(m):   # fraction of the mask's (non-DC) power beyond shell 20 (8 A)
    P = np.bincount(sh, (np.abs(np.fft.fftn(m)) ** 2).ravel(), n // 2 + 1); P[0] = 0; return P[20:].sum() / P.sum()
x = np.arange(n) - n // 2; rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); sph = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
lpc = lp(cons, 10.); env = (lpc > 0.1 * lpc.max()).astype(np.float32); envs = np.clip(lp(env, 10.), 0, 1)
d = A + '/fam1/zoomregionmotion_1_s05/level1/maps'; e, o = read_mrc(d + '/recvol_state01_even.mrc').astype(np.float32), read_mrc(d + '/recvol_state01_odd.mrc').astype(np.float32)
sel = lpc > 0.2 * lpc.max()
def dev(v): Aa = np.stack([lpc[sel], np.ones(sel.sum())], 1); ab = np.linalg.lstsq(Aa, v[sel], rcond=None)[0]; return (v - ab[1]) / ab[0] - lpc
de, do = dev(e), dev(o)
for name, fn in (('region_1', A + '/fam1/regions/region_1.mrc'), ('motion_1', A + '/fam1/regions_motion/region_1.mrc')):
    m = read_mrc(fn).astype(np.float32); print('%s: shape %s min %.3g max %.3g mean %.3g, fraction > 0.5: %.3f' % (name, m.shape, m.min(), m.max(), m.mean(), (m > 0.5).mean()))
    for tag, mk in (('raw', m), ('lp10', np.clip(lp(m, 10.), 0, 1)), ('lp10*envs', np.clip(lp(m, 10.), 0, 1) * envs), ('lp20', np.clip(lp(m, 20.), 0, 1)), ('lp20*envs', np.clip(lp(m, 20.), 0, 1) * envs)):
        c = fsc(de * mk, do * mk); print('  %-10s spectral power beyond 8 A: %.4f ; own FSC of the closed sub-state, shells 4..28 step 4: %s' % (tag, spec_frac(mk), ' '.join('%.2f' % v for v in c[4:29:4])))
for tag, mk in (('sphere', sph), ('envs', envs), ('env hard', env)):
    c = fsc(de * mk, do * mk); print('  %-10s spectral power beyond 8 A: %.4f ; own FSC: %s' % (tag, spec_frac(mk), ' '.join('%.2f' % v for v in c[4:29:4])))
c = fsc(e * sph, o * sph); print('  raw maps (not deviations), sphere: %s' % ' '.join('%.2f' % v for v in c[4:29:4]))
c = fsc(de, do); print('  no mask: %s' % ' '.join('%.2f' % v for v in c[4:29:4]))
print('DIAG MASK DONE')
PYEOF

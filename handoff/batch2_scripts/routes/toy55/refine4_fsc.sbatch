#!/bin/bash
#SBATCH -J klem_refine4_fsc
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# unfiltered half-map FSC of the refined states (recvol_stateXX_{even,odd}_unfil.mrc) vs the before maps
/mnt/nasapps/production/python/3.12/bin/python3 - <<'PYEOF'
import numpy as np, sys, os, glob
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc
A = '/mnt/beegfs/elmlund/afan/rnasep_klem/f692e'
def fsc_curve(fe, fo, smpd=1.464):
    e, o = read_mrc(fe).astype(np.float32), read_mrc(fo).astype(np.float32); n = e.shape[0]; x = np.arange(n) - n // 2
    rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
    E, O = np.fft.fftn(e * msk), np.fft.fftn(o * msk); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2)
    num = np.bincount(sh.ravel(), np.real(np.conj(E) * O).ravel()); d1 = np.bincount(sh.ravel(), np.abs(E).ravel() ** 2); d2 = np.bincount(sh.ravel(), np.abs(O).ravel() ** 2)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * smpd / k, fsc
rd = sorted(glob.glob(A + '/refine4/*_refine3D_states'))[-1]
before = {1: A + '/consolidated_12_fullw_region/fam0/regionall/level1/maps/recvol_state01', 2: A + '/consolidated_12_fullw_region/fam0/regionall/level1/maps/recvol_state16', 3: A + '/consolidated_12_fullw_region/fam0/regionall/level1/maps/recvol_state29', 4: A + '/consolidated_12_fullw_zoom/fam0/regionall/level1/maps/recvol_state01'}
print('state  before(A)  after_unfil(A)  after_filt(A)   [FSC at 12/10/8/6 A, unfiltered after]')
for s in range(1, 5):
    rb, _ = fsc_curve(before[s] + '_even.mrc', before[s] + '_odd.mrc')
    fu = '%s/recvol_state%02d_even_unfil.mrc' % (rd, s); ru, cu = fsc_curve(fu, fu.replace('_even', '_odd')) if os.path.isfile(fu) else (float('nan'), None)
    ff = '%s/recvol_state%02d_even.mrc' % (rd, s); rf, _ = fsc_curve(ff, ff.replace('_even', '_odd'))
    n = 160; sm = 1.464; vals = ' '.join('%.2f' % cu[int(round(n * sm / r))] for r in (12, 10, 8, 6)) if cu is not None else '-'
    print('%d  %.1f  %.1f  %.1f   [%s]' % (s, rb, ru, rf, vals))
print('FSC DONE')
PYEOF

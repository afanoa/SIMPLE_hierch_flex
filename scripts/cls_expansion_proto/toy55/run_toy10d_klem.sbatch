#!/bin/bash
#SBATCH -J toy10d_klem
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# toy10d after the simulation: refined consensus poses -> latd/ptcl3D_r0.txt, truth axis over the 11 states, then klem
# (random init, 50k) and a composition readout for conformer 63. env: NSUB (50000), KTAG (_fastrand), INIT (random)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; OUT=$W/latd; mkdir -p $OUT
[ -s $W/projd/ptcl3D_refined.txt ] || { echo "no refined pose table for toy10d"; exit 1; }
cp -f $W/projd/ptcl3D_refined.txt $OUT/ptcl3D_r0.txt
STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $W/rnasep_conformers10d.csv)
TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*_m160.mrc $W/maps160/variants_big/conf$(printf %02d $c)_*_m160.mrc 2>/dev/null | head -1; done | tr "\n" , | sed 's/,$//')
AXIS=$($PY - <<PYEOF
import numpy as np
def rd(fn):
    h = np.fromfile(fn, np.int32, 256); nx, ny, nz = h[:3]; return np.fromfile(fn, np.float32, offset=1024 + int(h[23]), count=nx*ny*nz).reshape(nz, ny, nx).astype(np.float64)
fns = '$TRUTH'.split(','); V = [rd(f) for f in fns]; B = V[0].shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx-B/2, yy-B/2), zz-B/2)
msk = r < 100/1.464; F = [np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(v * msk))) for v in V]; band = (np.rint(r) >= B*1.464/60) & (np.rint(r) <= B*1.464/10)
X = np.array([np.concatenate([f[band].real, f[band].imag]) for f in F]); X -= X.mean(0); U, s, Vt = np.linalg.svd(X, full_matrices=False); pc = U[:, 0] * s[0]; pc = pc / np.abs(pc).max()
print(','.join('%.3f' % v for v in pc))
PYEOF
); echo "truth axis (PC1 of the 11 states, STATES order $STATES): $AXIS"
export TAG=d AXIS=$AXIS INIT=${INIT:-random} NSUB=${NSUB:-50000} KTAG=${KTAG:-_fastrand} SIMDIR=simd PROJ=projd CSV=rnasep_conformers10d.csv
bash $W/run_klem.sbatch
D=$OUT/klem$KTAG; $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$D/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z1, z2, am, gm = L[:, 1], L[:, 2], L[:, 3].astype(int), L[:, 5]
import sys; sys.path.insert(0, '$W'); from latent_windows import read_gt
gt = read_gt('$W/simd/gt_particles.txt'); conf = np.array([gt.get(int(p), (0, 0))[0] for p in pind]); ncell = 49
def cell(z1v, z2v): return np.clip(np.rint(z1v), 0, 6).astype(int) * 7 + np.clip(np.rint(z2v), 0, 6).astype(int)
c = cell(z1, z2); on = gm > 0.5; m63 = conf == 63
print('toy10d composition readout: conformer 63 = %d particles (%.3f of conformer particles)' % (m63.sum(), m63.sum() / max((conf > 0).sum(), 1)))
print('  conformer 63 in the free (junk) components %.3f vs other conformers %.3f' % (np.mean(~on[m63]), np.mean(~on[(conf > 0) & ~m63])))
cnt = np.bincount(c[on & m63], minlength=ncell); tot = np.bincount(c[on & (conf > 0)], minlength=ncell); top = np.argsort(-cnt)[:5]
print('  cells holding conformer 63 (share of its on-grid particles / purity of the cell for 63): ' + ' '.join('(%d,%d): %.2f/%.2f' % (k // 7, k % 7, cnt[k] / max(cnt.sum(), 1), cnt[k] / max(tot[k], 1)) for k in top))
pur = cnt / np.maximum(tot, 1); print('  best cell purity for 63: %.2f (base rate %.2f); cells with purity > 0.5: %d' % (pur.max(), m63[on].mean(), int((pur > 0.5).sum())))
par = {int(l.split(',')[1]): int(l.split(',')[2]) for l in open('$W/rnasep_conformers10d.csv') if not l.startswith('emdb')}; parent = max(par, key=par.get)
mp = conf == parent; cp = np.bincount(c[on & mp], minlength=ncell)
print('  overlap with the parent conformer %d: fraction of 63 particles in cells where the parent is the majority %.2f' % (parent, float(cnt[(cp > cnt)].sum() / max(cnt.sum(), 1))))
PYEOF

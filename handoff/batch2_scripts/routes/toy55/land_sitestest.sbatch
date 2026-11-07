#!/bin/bash
#SBATCH -J land_sitestest
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/landscape/land_sitestest.%j.out
# discrete-state basin test on simsites: morph combination 2 (arm open) -> 6 (arm closed), pose-marginalised scores, basins vs truth
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; L=$K_/landscape; O=$L/sites_arm; mkdir -p $O; cd $O
log(){ echo "[$(date '+%F %T') sitestest] $*"; }
[ -s poses.txt ] || simple_exec prg=print_project_field projfile=$K_/rflex_simsites/node_r/fstageB/2_refine3D/b.simple oritype=ptcl3D mkdir=no > poses.txt 2> /dev/null; log "poses $(grep -c . poses.txt)"
for h in A B; do [ -s $K_/lsamp/caches/simsites_${h}_c48.npy ] || $PY $T/land_cache.py $T/simsites_$h $K_/lsamp/caches/simsites_${h}_c48.npy; done
$PY - <<PYEOF
import numpy as np
cmb = {int(l.split()[0]): int(l.split()[3]) for l in open('$T/maps160/simsites/combos.txt') if l.strip() and not l.startswith('#')}
g = np.loadtxt('$T/simsites/gt_particles.txt', comments='#'); conf = g[:, 1].astype(int)
np.savetxt('$O/rows_nonjunk.txt', (conf > 0).astype(int), fmt='%d')
with open('$O/gt_arm.txt', 'w') as f:
    f.write('# particle code(200 + 30 C_arm_closed; 0 junk) closed e1 e2\n')
    for i, c in enumerate(conf): f.write('%d %d 0 0 0\n' % (i + 1, 200 + 30 * cmb[c] if c > 0 else 0))
print('non-junk %d, arm closed %d' % ((conf > 0).sum(), sum(1 for c in conf if c > 0 and cmb[c] == 1)))
PYEOF
[ -s $O/morph/path_30.mrc ] || $PY -u $T/land_morph.py $T/maps160/simsites/conf02_sites.mrc $T/maps160/simsites/conf06_sites.mrc $O/morph 2>&1 | grep -v Warning
ls $O/morph/path_*.mrc > $O/maps.txt
[ -s $O/ps.npz ] || $PY -u $T/land_pscore.py --sim $T/simsites_A,$T/simsites_B --mr $T/simsites_A/meta_rows.txt,$T/simsites_B/meta_rows.txt --cache $K_/lsamp/caches/simsites \
   --poses $O/poses.txt --rowflags $O/rows_nonjunk.txt --maps $O/maps.txt --out $O/ps.npz --ctfabs 0 --mskdiam 232 2>&1 | grep -v "Warning\|block"
for lam in 0.001 0.01; do log "basins lam $lam"; $PY -u $T/land_basins.py --scores $O/ps.npz --gt $O/gt_arm.txt --shapes sim --lam $lam --thr 0.5,0.6,0.75 2>&1 | grep -v Warning; done
log DONE

#!/bin/bash
#SBATCH -J rnasep_pairing_check
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Does the re-imported project pair the same image with the same pose as the meta/klem path? For 8 particles of the
# within-10k set: the image klem reads (stkind/indstk of the meta) vs the image the imported project row points to.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; P=$A/rnasep_klem/expand1/fam1_pipeline; KD=$P/klem_fb; META=$P/meta_fb.npz; IC=$KD/region_align_10k/import_check
simple_exec prg=print_project_field projfile=$IC/plain/proj/proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $IC/plain/ptcl3D_imported.txt
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; simple_exec prg=print_project_field projfile=$PS oritype=ptcl3D mkdir=no 2>/dev/null | head -400000 > $IC/family_ptcl3D_head.txt
$PY - <<PYEOF
import numpy as np, glob
M = np.load('$META'); stk = M['stkind'].astype(int); ind = M['indstk'].astype(int); e1 = M['e1']; N = len(stk); order = np.lexsort((ind, stk)); pos = np.empty(N, int); pos[order] = np.arange(N)
files = sorted(glob.glob('/mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs'))
def img(f, i):
    h = np.fromfile(f, np.int32, 256); nx, ny = int(h[0]), int(h[1]); off = 1024 + int(h[23]) + (i - 1) * nx * ny * 4; return np.fromfile(f, np.float32, nx * ny, offset=off)
rows = [dict(t.split('=', 1) for t in l.split() if '=' in t) for l in open('$IC/plain/ptcl3D_imported.txt') if l.strip()]
fam = [dict(t.split('=', 1) for t in l.split() if '=' in t) for l in open('$IC/family_ptcl3D_head.txt') if l.strip()]
sub = np.loadtxt('$KD/closed_focal_docked30k/rank/subset_top10k.txt', dtype=int)[:8]
print('imported project rows: %d; keys of row 1: stkind %s indstk %s; family project head rows: %d' % (len(rows), rows[0].get('stkind'), rows[0].get('indstk'), len(fam)))
for r in sub:
    i = r - 1; a = img(files[stk[i] - 1], ind[i]); m = pos[i]; rr = rows[m]; sk, ik = int(float(rr['stkind'])), int(float(rr['indstk'])); b = img(files[sk - 1], ik)
    famrow = fam[i] if i < len(fam) else {}
    print('row %d: meta (stk %d, ind %d) e1 %.2f | imported row %d (stk %d, ind %d) e1 %s state %s | family project row e1 %s stk %s ind %s | image corr %.3f' % (r, stk[i], ind[i], e1[i], m + 1, sk, ik, rr.get('e1'), rr.get('state'), famrow.get('e1'), famrow.get('stkind'), famrow.get('indstk'), float(np.corrcoef(a, b)[0, 1])))
PYEOF

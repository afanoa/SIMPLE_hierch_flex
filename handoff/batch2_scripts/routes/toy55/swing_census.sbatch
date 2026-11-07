#!/bin/bash
#SBATCH -J klem_swing_census
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 10:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# census on swing30 subsets (S1 continuum 0-30 deg, S2 three clusters 0-2 / 14-16 / 28-30, S3 14-16 + 8 % at 28-30): truth, swing mask, mean, focused flex per subset, census.py; env SETS OUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-csbdevel} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; S=$T/simswing30; P=$T/projswing30/projswing30.simple
O=${OUT:-$K_/factorial/census_swing}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') swing census] $*"; }
[ -s $O/poses.txt ] || simple_exec prg=print_project_field projfile=$P oritype=ptcl3D mkdir=no 2>/dev/null > $O/poses.txt
$PY - <<PYEOF || { log "setup FAILED"; exit 1; }
import numpy as np, os
G = np.loadtxt('$S/gt_particles.txt', comments='#', usecols=(0, 1)); N = int(G[:, 0].max()); conf = np.zeros(N, int); conf[G[:, 0].astype(int) - 1] = G[:, 1].astype(int)
sw = (conf >= 200) & (conf <= 230); ang = np.where(sw, conf - 200, -1); rng = np.random.default_rng(11)
sets = {'S1': sw, 'S2': np.isin(conf, [200, 201, 202, 214, 215, 216, 228, 229, 230]), 'S3': np.isin(conf, [214, 215, 216])}
rare = np.flatnonzero(np.isin(conf, [228, 229, 230])); nr = int(round(0.08 / 0.92 * sets['S3'].sum())); sets['S3'][rng.choice(rare, nr, replace=False)] = True
with open('$O/truth.txt', 'w') as f:
    f.write('# row angle junk\n')
    for i in range(N): f.write('%d %s %d\n' % (i + 1, ('%d' % ang[i]) if sw[i] else 'nan', 0 if sw[i] else 1))
def rd(p):
    h = np.fromfile(p, np.int32, 256); n = int(h[0]); return np.fromfile(p, np.float32, offset=1024 + int(h[23]), count=n ** 3).reshape(n, n, n)
def wr(p, x, ref):
    with open(ref, 'rb') as f: h = f.read(1024)
    hh = np.frombuffer(h, np.int32).copy(); hh[23] = 0; hf = hh.view(np.float32); hf[19:22] = x.min(), x.max(), x.mean()
    with open(p, 'wb') as f: f.write(hh.tobytes()); f.write(x.astype(np.float32).tobytes())
mp = {c: '$S/maps/conf%d_swing_m160.mrc' % c for c in range(200, 231)}
for k, v in sets.items():
    np.savetxt('$O/rows_%s.txt' % k, v.astype(int), fmt='%d'); cs, cn = np.unique(conf[v], return_counts=True)
    wr('$O/mean_%s.mrc' % k, sum(n_ * rd(mp[c]).astype(np.float64) for c, n_ in zip(cs, cn)) / cn.sum(), mp[200])
    print('%s: %d particles, angles %s' % (k, v.sum(), ' '.join('%d:%d' % (c - 200, n_) for c, n_ in zip(cs, cn))))
open('$O/swing_list.txt', 'w').write(''.join('c%d %s %s 1000\n' % (c, mp[c], mp[c]) for c in range(200, 231)))
PYEOF
[ -s $O/swing_mask.mrc ] || $PY $T/arm_mask2.py $O/swing_mask.mrc --list $O/swing_list.txt --mskdiam 232 --frac 0.2 --dil 3 --falloff 6 2>&1 | grep -v Warning | tail -1
for k in ${SETS:-S1 S2 S3}; do F=$O/flex_$k; mkdir -p $F; cd $F
  if [ -z "$(ls $F/*_flex_pca/flex_pca_coordinates.txt 2>/dev/null)" ]; then cp -f $P f.simple && simple_exec prg=selection projfile=f.simple oritype=ptcl3D infile=$O/rows_$k.txt mkdir=no > sel.log 2>&1 \
    && simple_exec prg=update_project projfile=f.simple qsys_name=slurm nparts=8 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=flex_pca projfile=f.simple vol1=$O/mean_$k.mrc mskdiam=232 lp=15 npreimages=8 rec_states=no min_state_frac=0 umap=no n_probe_iters=2 pcg_mskfile=$O/swing_mask.mrc nparts=8 nthr=8 > flex.log 2>&1 \
    || { log "$k flex FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-150)"; continue; }; fi
  CO=$(ls -td $F/*_flex_pca | head -1)/flex_pca_coordinates.txt; log "$k: focused flex done -> census"; cd $O
  $PY -u $T/census.py --full 1 --sim $S --poses $O/poses.txt --cache $T/benchcache_full_swing30 --ctfabs 0 --rows $O/rows_$k.txt --mask $O/swing_mask.mrc --mean $O/mean_$k.mrc \
    --coords $CO --truth $O/truth.txt --tcols angle --out $O/census_$k --pars ${PARS:-4} 2>&1 | grep -v Warning | tee $O/census_$k.txt; done
log DONE

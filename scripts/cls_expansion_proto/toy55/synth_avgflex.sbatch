#!/bin/bash
#SBATCH -J synth_avgflex
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Stitching by flex_pca on the direction-bin sub-class averages (env MODE; averages + directions from synth_bins_route):
# the averages imported as a particle set with their bin directions (class = bin), flex_pca in cavgs mode with VIEW
# CENTRING (the between-view part removed before the fit), GMM placement into K states, sub-class -> particle labels,
# truth readout (per-arm std ratio / per-site occupancy). Even/odd averages are fitted too as the cross-half gate.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_master1013/snap/c6e344c6b_reseed27; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16   # the flex build with flex_cavgs + flex_view_centre (pv3pm)
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; K=${K:-6}; NE=${NE:-10}; LP=${LP:-12}
D=$W/binsroute_$MODE; S=$D/stitch; O=$D/avgflex_k$K; rm -rf $O; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') avgflex $MODE] $*"; }
GT=$W/sim$MODE/gt_particles.txt; NPART=$(grep -vc '^#' $GT); CONS=$W/maps160/synth_$MODE/consensus_synth.mrc
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); M = np.load('$S/avgs_meta.npz'); N = int(M['N']); e1, e2, par = M['e1'], M['e2'], M['parent'].astype(int)
h = np.frombuffer(open('$S/avgs.mrcs', 'rb').read(1024), np.int32); nsym = int(h[23]); nx = int(h[0]); A = np.memmap('$S/avgs.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(3 * N, nx, nx))
for name, sl in (('full', slice(0, N)), ('even', slice(N, 2 * N)), ('odd', slice(2 * N, 3 * N))):
    ii = h.copy(); ii[2] = N; ii[9] = N; ii[23] = 0; X = np.ascontiguousarray(A[sl]); ff = ii.view(np.float32); ff[19] = float(X.min()); ff[20] = float(X.max()); ff[21] = float(X.mean())
    open('$O/avgs_%s.mrcs' % name, 'wb').write(ii.tobytes() + X.tobytes())
with open('$O/avgs_oritab.txt', 'w') as f:
    for g in range(N): f.write('e1=%.4f e2=%.4f e3=0.0000 x=0.0000 y=0.0000 class=%d state=1\n' % (e1[g], e2[g], par[g]))
print('%d averages x 3 stacks written, class = bin (%d bins)' % (N, len(set(par))))
PYEOF
for h in full even odd; do
  P=$O/run_$h; mkdir -p $P; cd $P; simple_exec prg=new_project projname=avg dir=$P mkdir=no > new_project.log 2>&1 || { log "new_project FAILED ($h)"; exit 1; }
  simple_exec prg=import_particles projfile=avg.simple stk=$O/avgs_$h.mrcs oritab=$O/avgs_oritab.txt ctf=no smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED ($h): $(grep -a -m2 ERROR import.log | cut -c1-160)"; exit 1; }
  log "$h: flex_pca on the averages: cavgs model, view-centred, lp $LP neigs $NE, gmm placement into $K states"
  simple_exec prg=flex_pca projfile=avg.simple vol1=$CONS mskdiam=230 lp=$LP neigs=$NE n_probe_iters=2 npreimages=$K state_placement=gmm flex_cavgs=model flex_view_centre=yes min_neff=5 rec_states=yes nparts=1 nthr=16 mkdir=no > flex.log 2>&1 || { log "flex_pca FAILED ($h): $(grep -a -m3 -i 'error\|THROW' flex.log | cut -c1-200)"; continue; }
  FD=$(dirname $(find $P -name flex_pca_coordinates.txt | head -1)); [ -n "$FD" ] && [ -s $FD/flex_pca_coordinates.txt ] || { log "no flex_pca_coordinates.txt under $P"; continue; }; grep -a -i -E "rows used|pruned|neff|state populations|placement|effective rank" $P/flex.log | head -6 | cut -c1-200
  awk '!/^#/ {print $1, $3}' $FD/flex_pca_coordinates.txt > $O/sublab_$h.txt; log "$h: $(awk '{print $2}' $O/sublab_$h.txt | sort | uniq -c | tr '\n' ' ')"
done
cd $O; $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$D/p/cls_expansion_labels.txt', comments='#').astype(int); sub = L[:, 2]
for h in ('full', 'even', 'odd'):
    try: sl = np.loadtxt('$O/sublab_%s.txt' % h, dtype=int)
    except Exception: continue
    lab = np.zeros(sl[:, 0].max() + 1, int); lab[sl[:, 0]] = sl[:, 1]; pl = np.zeros($NPART, int); pl[L[:, 0] - 1] = lab[np.minimum(sub, len(lab) - 1)]; np.savetxt('$O/particle_states_%s.txt' % h, pl, fmt='%d')
# cross-half agreement of the sub-class states (even vs odd fits, best label permutation)
try:
    e, o = np.loadtxt('$O/sublab_even.txt', dtype=int), np.loadtxt('$O/sublab_odd.txt', dtype=int); from itertools import permutations
    K = int(max(e[:, 1].max(), o[:, 1].max())); best = max(np.mean(np.array([p[x - 1] for x in e[:, 1]]) == o[:, 1]) for p in permutations(range(1, K + 1))) if K <= 7 else float('nan')
    print('even/odd fits: sub-class state agreement %.2f (best permutation of %d labels; chance %.2f)' % (best, K, 1. / K))
except Exception as ex: print('even/odd comparison skipped:', ex)
PYEOF
for h in full even odd; do [ -s $O/particle_states_$h.txt ] || continue; echo "== truth readout, $h averages"; $PY $W/label_truth_readout.py $O/particle_states_$h.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -$((K + 2)); done
log "DONE"

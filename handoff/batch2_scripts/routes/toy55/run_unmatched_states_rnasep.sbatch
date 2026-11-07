#!/bin/bash
#SBATCH -J rnasep_unmatched_states
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The reproducible real states that match none of the ten deposited conformers, as maps: leaves = family 1's two
# level-3 grandchildren + family 2's feedback board (no double counting), clustered at similarity SIM; clusters whose
# dominant deposited match is below 0.3 and that hold >= NMIN particles are reconstructed at full box at their family's
# refined poses (one reconstruct3D per family), scored against the ten deposited maps (set mean subtracted) and FSC.
# env: SIM (0.6), NMIN (5000), MAXK (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
D=$A/rnasep_klem/unmatched_states; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') un] $*"; }
LEAVES="$E/fam1_l3/sub1 $E/fam1_l3/sub2 $E/fam2_pipeline/klem_fb"; for d in $E/fam2_l3/sub1 $E/fam2_l3/sub2; do [ -s $d/klem_halves.npz ] && LEAVES="$LEAVES $d"; done
for d in $LEAVES; do [ -s $d/cell_matches.txt ] || $PY $W/deposited_match.py $d/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $d/cell_matches.txt > /dev/null 2>&1; done
$PY $W/count_states.py $LEAVES --save $D/clusters.txt --at ${SIM:-0.6} 2>&1 | grep -v Warning
$PY - <<PYEOF
import numpy as np, collections
C = [l.split() for l in open('$D/clusters.txt') if not l.startswith('#')]; N = len(np.load('$META')['flags']); fam_of = lambda d: 2 if '/fam2' in d else 1
clus = collections.defaultdict(list)
for b, d, cell, k, w, dep, corr in C: clus[int(k)].append((d, int(cell), float(w), int(dep), float(corr)))
rows = []
for k, cells in clus.items():
    wsum = sum(c[2] for c in cells); best = max(cells, key=lambda c: c[4]); fam = fam_of(cells[0][0])
    rows.append((k, wsum, best[3], best[4], fam, cells))
unm = sorted([r for r in rows if r[3] < 0.3 and r[1] >= ${NMIN:-5000}], key=lambda r: -r[1])[:${MAXK:-10}]
print('clusters at similarity ${SIM:-0.6}: %d; unmatched (best deposited corr < 0.3) with >= ${NMIN:-5000} particles: %d' % (len(rows), len(unm)))
for fam in (1, 2):
    st = np.zeros(N, int); names = []
    for j, r in enumerate([r for r in unm if r[4] == fam]):
        for d, cell, w, dep, corr in r[5]:
            L = np.loadtxt(d + '/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int) - 1; gm = L[:, 5]; st[pind[(gm > 0.5) & (am == cell)] - 1] = j + 1
        names.append('cluster %d (%d cells, %.0f particles, best deposited %d at %+.2f)' % (r[0], len(r[5]), r[1], r[2], r[3]))
    np.savetxt('$D/states_fam%d.txt' % fam, st, fmt='%d'); open('$D/names_fam%d.txt' % fam, 'w').write('\n'.join(names) + '\n'); print('family %d: %d unmatched states, %d particles: %s' % (fam, len(names), int((st > 0).sum()), ' | '.join(names)))
PYEOF
for fam in 1 2; do
  K=$(grep -c . $D/names_fam$fam.txt); [ "$K" -eq 0 ] && continue; P=$E/fam${fam}_pipeline; PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; mkdir -p $D/fam$fam && cd $D/fam$fam && cp -f $PS p.simple
  simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states_fam$fam.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED (family $fam)"; cd $D; continue; }
  log "family $fam: reconstruct3D of $K unmatched states at the refined poses"
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$K nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; cd $D; continue; }
  RD=$(ls -td $D/fam$fam/*_reconstruct3D | head -1); mkdir -p $D/fam$fam/maps; for w in $(seq -f %02g 1 $K); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/fam$fam/maps/ 2>/dev/null; done; done
  $PY $W/score_maps.py $D/fam$fam/maps --n $K 2>&1 | grep -A $K 'conformer - consensus, 60-20' | head -$((K + 1)) | sed "s/^/[fam $fam] /"
  $PY - <<PYEOF
import numpy as np, sys, glob
sys.path.insert(0, '$W'); from model_class import read_mrc, Shells, fftc
S = Shells(160, 1.464); sh = S.sh
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return 160 * 1.464 / s
names = [l.strip() for l in open('$D/names_fam$fam.txt') if l.strip()]
for k, nm in enumerate(names):
    try: print('  [fam $fam] map %2d: FSC 0.143 %.1f A | %s' % (k + 1, fsc_res(read_mrc('$D/fam$fam/maps/recvol_state%02d_even.mrc' % (k + 1)), read_mrc('$D/fam$fam/maps/recvol_state%02d_odd.mrc' % (k + 1))), nm))
    except Exception as e: print('  [fam $fam] map %2d: %s (%s)' % (k + 1, nm, e))
PYEOF
  cd $D
done
log DONE

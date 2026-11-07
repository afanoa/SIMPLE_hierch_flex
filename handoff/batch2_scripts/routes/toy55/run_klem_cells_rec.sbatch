#!/bin/bash
#SBATCH -J rnasep_klem_rec
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Full-box reconstructions of bins along the main axis of a klem latent on real RNase P (SIMPLE reconstruct3D), then a
# comparison of every bin map with the 10 docked Degenhardt conformers (same frame as the toy). env: KDIR (klem dir), NB (7), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KDIR=${KDIR:-$A/rnasep_klem/klem_full_npz}; NB=${NB:-7}; NP=${NP:-10}
SRC=$A/rnasep_coupled/full692k_2d/rnasep.simple; META=$A/rnasep_coupled/basis_span_all/meta_all.npz; D=$KDIR/rec${RTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') rec] $*"; }
if [ -n "${STATES:-}" ]; then cp -f $STATES $D/states.txt; log "states from $STATES: $(sort -n $D/states.txt | uniq -c | awk '$2>0' | tr '\n' ' ')"; else $PY $W/klem_bins.py $KDIR $META $D/states.txt --nb $NB --pc ${PC:-1}; fi
cp -f $SRC p.simple
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; tail -2 sel.log; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=${MSKDIAM:-200} nstates=$NB nparts=$NP nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in $(seq -f %02g 1 $NB); do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
log "maps: $(ls $D/maps | wc -l) files; cFAR $(grep -a cFAR rec.log | awk '{print $NF}' | tr '\n' ' ')"
# compare with the docked Degenhardt conformers (toy frame = consensus frame): corr(bin - mean of bins, conformer - toy consensus) per band
STATES=1,10,18,27,28,37,46,48,49,51; TRUTH=$(for c in $(echo $STATES | tr , " "); do ls $W/maps160/docked/conf$(printf %02d $c)_*.mrc | head -1; done | tr "\n" , | sed 's/,$//')
awk -v nb=$NB 'BEGIN{for(i=1;i<=nb;i++) print i, 0, 0, 0, 0, 0, 0}' > $D/windows.txt
$PY - $D $NB "$TRUTH" "$STATES" <<'PYEOF'
import sys, os, numpy as np
sys.path.insert(0, '/mnt/beegfs/elmlund/afan/routes/toy55'); from model_class import read_mrc, Shells, fftc, BANDS
D, NB, truth, states = sys.argv[1], int(sys.argv[2]), sys.argv[3].split(','), [int(s) for s in sys.argv[4].split(',')]
S = Shells(160, 1.464); msk = S.soft_mask(100 / 1.464); cons = read_mrc('/mnt/beegfs/elmlund/afan/routes/toy55/maps160/toy_consensus.mrc') * msk
T = {s: fftc(read_mrc(f) * msk - cons) for s, f in zip(states, truth)}
M = [read_mrc('%s/maps/recvol_state%02d.mrc' % (D, w)) * msk for w in range(1, NB + 1)]; Mm = np.mean(M, 0); sh = S.sh
def bc(x, y, lo, hi):
    m = S.band(lo, hi); return float(np.real(np.sum(np.conj(x[m]) * y[m])) / max(np.sqrt(np.sum(np.abs(x[m]) ** 2) * np.sum(np.abs(y[m]) ** 2)), 1e-30))
def fsc_res(x, y):
    a_, b_ = fftc(x), fftc(y); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return 160 * 1.464 / s
print('bin maps vs the 10 docked Degenhardt conformers: corr(bin - mean, conformer - consensus) at 60-20 A; best conformer per bin')
for w in range(NB):
    Dw = fftc(M[w] - Mm); cc = {s: bc(Dw, T[s], 60, 20) for s in states}; best = max(cc, key=cc.get)
    try: res = fsc_res(read_mrc('%s/maps/recvol_state%02d_even.mrc' % (D, w + 1)) * msk, read_mrc('%s/maps/recvol_state%02d_odd.mrc' % (D, w + 1)) * msk)
    except Exception: res = 0.
    print('  bin %d: FSC 0.143 %.1f A | best conformer %d (%.2f) | all: %s' % (w + 1, res, best, cc[best], ' '.join('%d:%+.2f' % (s, cc[s]) for s in states)))
# ends: corr of (bin NB - bin 1) with each conformer pair difference
E = fftc(M[-1] - M[0]); print('end-to-end difference vs conformer differences at 60-20 A: ' + ' '.join('%d-%d:%+.2f' % (a, b, bc(E, T[a] - T[b], 60, 20)) for a in (27, 28, 49) for b in (10, 37, 18)))
PYEOF
log DONE

#!/bin/bash
#SBATCH -J rnasep_closed_zoom
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# A dedicated closed-state board: the particles of family 1's feedback board whose argmax cell matches a closed-family
# conformer (18, 46, 10), at the refined poses, focus mask on (the reproducible variability inside the closed family =
# the arm), pose marginalisation, balanced E-step, 12->8 A, 7x7; cells scored by closed-arm occupancy and by match to
# 18; the top cells' particles reconstructed at full box at the refined poses -> occupancy / FSC of the closed map.
# env: NTOP (3), BAL (0.5), KD (board), CONFS (18,46,10), MSKD (232: reference mask diameter, the default 200 clips the
#      arm tip at 100 A radius), LAYERS ("--nlayer 2 --tie_layers" = occupancy model: one shared additive component), ZTAG
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KD=${KD:-$P/klem_fb}; META=$P/meta_fb.npz; METAALL=$A/rnasep_coupled/basis_span_all/meta_all.npz
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
D=$KD/closed_zoom${ZTAG:-}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') cz] $*"; }
$PY - <<PYEOF
import numpy as np
M = np.loadtxt('$KD/matched_rec/cell_matches.txt', comments='#'); best = {int(r[0]): int(r[1]) for r in M}; confs = [int(c) for c in '${CONFS:-18,46,10}'.split(',')]
L = np.loadtxt('$KD/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]; sel = (gm > 0.5) & np.array([best.get(int(c), 0) in confs for c in am])
np.savetxt('$D/closed_subset.txt', pind[sel], fmt='%d'); print('closed-family subset: %d particles from %d cells matched to %s' % (sel.sum(), sum(1 for v in best.values() if v in confs), confs))
PYEOF
export OMP_NUM_THREADS=1; log "closed board at the refined poses"
$PY $W/klem.py $STK none none $D/board --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $D/closed_subset.txt --ctf auto \
    --lps "12 12 10 10 8 8 8 8 8" --sigmas "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" --mskdiam ${MSKD:-232} ${LAYERS:-} --focus 1 --pose_marg 3 --devw 0.05 --pnoise --balance ${BAL:-0.5} --tag "rnasep closed zoom${ZTAG:-}" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "subset\|focus mask\|final:\|states: cross-half\|Error\|Traceback"
echo "== identification"; $PY $W/deposited_match.py $D/board/cells $W/maps160/docked $W/maps160/toy_consensus.mrc --save $D/board/cell_matches.txt | tail -13
echo "== closed-arm occupancy per cell"; $PY $W/closed_arm_cells.py $D/board --top 8
$PY - <<PYEOF
import numpy as np
C = np.loadtxt('$D/board/closed_arm_cells.txt', comments='#'); order = np.argsort(-C[:, 1]); top = [int(C[k, 0]) for k in order[:${NTOP:-3}]]
L = np.loadtxt('$D/board/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); gm = L[:, 5]; N = len(np.load('$METAALL')['flags']); st = np.zeros(N, int)
for j, c in enumerate(top): st[pind[(gm > 0.5) & (am == c)] - 1] = 1   # one state: the union of the top cells
st2 = st.copy(); st2[pind[(gm > 0.5) & (am == top[0])] - 1] = 2         # state 2 = the single best cell (state 1 keeps the rest of the top cells)
np.savetxt('$D/states.txt', st2, fmt='%d'); print('top cells by closed-arm occupancy: %s -> state 1 (union minus best) %d particles, state 2 (best cell) %d particles' % (top, int((st2 == 1).sum()), int((st2 == 2).sum())))
PYEOF
export OMP_NUM_THREADS=8; PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple; cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$D/states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=240 nstates=2 nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $D/*_reconstruct3D | head -1); mkdir -p $D/maps; for w in 01 02; do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc $D/maps/ 2>/dev/null; done; done
$PY - <<PYEOF
import numpy as np, glob, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, Shells, fftc
from scipy.ndimage import gaussian_filter
W = '$W'; SMPD = 1.464; S = Shells(160, SMPD); sh = S.sh; msk = S.soft_mask(110 / SMPD); g = lambda v: gaussian_filter(v, 2.5)
cons = read_mrc(W + '/maps160/toy_consensus.mrc'); c18 = read_mrc(glob.glob(W + '/maps160/docked/conf18_*.mrc')[0]); L18, Lc = g(c18), g(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); diff = L18 - Lc
arm = (diff > diff.mean() + 3 * diff.std()) & dens18; core = dens18 & (np.abs(diff) < 0.5 * diff.std()); ref = L18[arm].mean() / L18[core].mean()
def fsc_res(x, y):
    a_, b_ = fftc(x * msk), fftc(y * msk); num = np.bincount(sh.ravel(), np.real(np.conj(a_) * b_).ravel()); d1 = np.bincount(sh.ravel(), np.abs(a_.ravel()) ** 2); d2 = np.bincount(sh.ravel(), np.abs(b_.ravel()) ** 2)
    f = num / np.maximum(np.sqrt(d1 * d2), 1e-30); s = 1
    while s < 80 and f[s] > 0.143: s += 1
    return 160 * SMPD / s
st = np.loadtxt('$D/states.txt', dtype=int)
for w, name in ((1, 'top cells minus the best'), (2, 'best cell')):
    v = read_mrc('$D/maps/recvol_state%02d.mrc' % w); occ = (g(v)[arm].mean() / g(v)[core].mean()) / ref
    try: res = fsc_res(read_mrc('$D/maps/recvol_state%02d_even.mrc' % w), read_mrc('$D/maps/recvol_state%02d_odd.mrc' % w))
    except Exception: res = float('nan')
    print('CLOSED ZOOM map %d (%s, %d particles): closed-arm occupancy %.2f | FSC 0.143 %.1f A' % (w, name, int((st == w).sum()), occ, res))
PYEOF
log DONE

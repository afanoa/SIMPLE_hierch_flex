#!/bin/bash
#SBATCH -J synth_cellllr
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# whitened per-image likelihood (arm_llr.py, fitted scale, true poses) of the images against a board's OWN cell maps, truth =
# the oracle start label: does the whitened scorer put the particles back in their start cells where the board's E-step does
# not? env BOARD (oracle dir with init_positions/init_labels and board/cells), MODE, REG, LP (20), NSUB (3000)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; BOARD=${BOARD:?}; MODE=${MODE:?}; REG=${REG:?}
O=$BOARD/cell_llr_lp${LP:-20}${HALF:+_$HALF}${HP:+_hp$HP}${LSSCALE:+_ls$LSSCALE}; mkdir -p $O/sim $O/maps; export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
ln -sf $W/sim$MODE/all.txt $O/sim/all.txt; ln -sf $W/sim$MODE/all.mrcs $O/sim/all.mrcs
$PY - <<PYEOF
import numpy as np, os, collections
P = np.loadtxt('$BOARD/init_positions.txt', dtype=int); L = np.loadtxt('$BOARD/init_labels.txt', dtype=int)
cell = {}; votes = collections.defaultdict(collections.Counter)
for pind, g1, g2 in P:
    if L[pind - 1] > 0: votes[L[pind - 1]][(g1, g2)] += 1
for lab, c in votes.items(): cell[lab] = c.most_common(1)[0][0]
with open('$O/sim/gt_particles.txt', 'w') as f:
    f.write('# pind label (oracle start cell)\n')
    for pind, g1, g2 in P: f.write('%d %d\n' % (pind, L[pind - 1]))
for lab, (g1, g2) in cell.items():
    src = '$BOARD/board/cells/cell_%02d_%02d.mrc' % (g1, g2); dst = '$O/maps/conf%02d_synth.mrc' % lab
    if '${HALF:-}':   # score against the half-set cell (cells_half/cell_NN_<half>.mrc): NN = the half cell most correlated with this merged cell
        import glob
        rd = lambda f: np.fromfile(f, np.float32, offset=1024)
        m = rd(src); best = max(glob.glob('$BOARD/board/cells_half/cell_*_${HALF:-}.mrc'), key=lambda f: float(np.corrcoef(m, rd(f))[0, 1])); src = best
    if os.path.exists(dst): os.remove(dst)
    os.symlink(src, dst)
print('%d labels -> cells, %d particles' % (len(cell), len(P)))
PYEOF
echo "== $(basename $BOARD): whitened likelihood of the images against the board's cells, truth = start cell"
$PY $W/arm_llr.py $O/sim none $REG $O/maps $O --nsub ${NSUB:-3000} --lp ${LP:-20} ${HP:+--hp $HP} ${LSSCALE:+--lsscale $LSSCALE} 2>&1 | grep -v Warning | grep -a -v "scored [0-9]"
echo "[$(date '+%F %T') cellllr] DONE"

#!/bin/bash
#SBATCH -J lrec_board
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/lrec30/logs/%x.%j.out
# converged-board test on the oracle swing family (simswing30, realistic SNR): klem boards (Fortran, local, near-true poses,
# random start) with the fast-loop schedule (5 rounds) and a converged one (15 rounds, annealed lattice sigma), on a 1-D chain
# (N1 x 1) and the usual 7 x 7; readout = the posterior-mean position (z1) and the argmax cell against the true swing.
# env LABEL (board), N1 (16), N2 (1), LPS, SIGMAS
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; KD=$A/rnasep_klem
export FB=$KD/f692e/klem2_bin CTFABS=0 PSCALE=0 SERVER_LIST= NPARTS= NTHR=16 INIT=${INIT:-random} NFREE=0
O=$KD/lrec30/board_${LABEL:?}; mkdir -p $O && cd $O; GT=$T/simswing30/gt_particles.txt
awk '$1 == 1 {print NR}' $KD/ff2_swing30/oracle_swing_rows.txt > fam_subset.txt
cp -f $T/projswing30/projswing30.simple p.simple
N1=${N1:-16} N2=${N2:-1} bash $T/fboard.sh $O/board $O/p.simple $O/fam_subset.txt "${LPS:?}" "${SIGMAS:?}" > board.out 2>&1 || { echo "board FAILED: $(tail -3 board.out)"; exit 1; }
tail -3 board.out
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); z1, z2, am, mr = L[:, 1], L[:, 2], L[:, 3].astype(int), L[:, 4]
gt = np.loadtxt('$GT', comments='#', usecols=(1,), dtype=int)[pind - 1]; x = gt - 200.; sw = (gt >= 200) & (gt <= 230)
def eta2(lab):
    W = sum(((x[sw & (lab == v)] - x[sw & (lab == v)].mean()) ** 2).sum() for v in np.unique(lab[sw])); return 1 - W / ((x[sw] - x[sw].mean()) ** 2).sum()
r1, r2 = np.corrcoef(z1[sw], x[sw])[0, 1], np.corrcoef(z2[sw], x[sw])[0, 1]; c = z1 if abs(r1) >= abs(r2) else z2; s = np.sign(np.corrcoef(c[sw], x[sw])[0, 1])
q = np.quantile(c, np.linspace(0, 1, 9)[1:-1]); b = np.searchsorted(q, c); o = np.argsort(s * c); n = len(o)
print('BOARD $LABEL: %d particles, median maxr %.3f; corr(z1, swing) %.3f, corr(z2, swing) %.3f; argmax cells eta2 %.3f; 8 quantile bins of the best z: eta2 %.3f, bin means %s; bottom/top 2%% %.1f/%.1f, 5%% %.1f/%.1f'
      % (len(pind), np.median(mr), abs(r1), abs(r2), eta2(am), eta2(b), ' '.join('%.1f' % x[b == v].mean() for v in (range(8) if s > 0 else range(7, -1, -1))), x[o[:n // 50]].mean(), x[o[-n // 50:]].mean(), x[o[:n // 20]].mean(), x[o[-n // 20:]].mean()))
PYEOF
echo DONE

#!/bin/bash
#SBATCH -J klem_gridding_ab
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# A/B of Kaiser-Bessel gridding + deapodisation in the board (env GRIDDING 0|1, MODE synthetic set, NSUB particles):
# one parent-level board at the refined poses on a random subset, same seed; readouts = cross-half reproducibility of
# the cells, half-map FSC of the 5 most populated cells, truth readout of the argmax cells (per-arm / per-site).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:-tethered2}; G=${GRIDDING:-0}; R=$A/rnasep_klem/synth_$MODE; O=$A/rnasep_klem/gridding_ab/${MODE}_g$G; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') gridab $MODE g$G] $*"; }
FB=${FB:-$A/rnasep_klem/gridding_ab/klem2_bin}; export FB SIMPLE_QSYS=local
GT=$W/sim$MODE/gt_particles.txt; NP=$(grep -vc '^#' $GT)
$PY -c "import numpy as np; rng=np.random.default_rng(0); n=$NP; s=np.sort(rng.choice(n, ${NSUB:-50000}, replace=False)+1); np.savetxt('$O/subset.txt', s, fmt='%d'); print('subset', len(s))"
log "board on $(grep -c . $O/subset.txt) particles, gridding $G ($(basename $FB))"
SERVER_LIST= NPARTS= NTHR=16 N1=7 N2=7 NFREE=1 CTFABS=1 PSCALE=1 INIT=random WRITE_HALVES=1 GRIDDING=$G bash $W/fboard.sh $O/board $R/src.simple $O/subset.txt "20 14 12 10 10 8 8 8" "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" > $O/board.out 2>&1 || { log "board FAILED: $(tail -3 $O/board.out | cut -c1-160)"; exit 1; }
grep -a "gridding\|Kaiser\|timing" $O/board.out | tail -3 | cut -c1-160
echo "== cell reproducibility (klem_states.txt: cell population reproducibility)"; awk '!/^#/ {n++; s+=$3; if($3>0.9) h++} END {printf "%d cells, mean cross-half reproducibility %.3f, %d cells > 0.9\n", n, s/n, h}' $O/board/klem_states.txt
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W')
H = np.load('$O/board/klem_halves.npz'); keys = list(H.keys()); print('halves keys:', keys[:6])
PYEOF
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); lab = np.zeros($NP, int); lab[pind - 1] = am; np.savetxt('$O/cells.txt', lab, fmt='%d')
PYEOF
echo "== truth readout of the argmax cells"; $PY $W/label_truth_readout.py $O/cells.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -6
log DONE

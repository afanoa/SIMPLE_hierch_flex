#!/bin/bash
#SBATCH -J synth_oracle
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# oracle test of the focused region board on a hinged-arm synthetic set (env MODE, REG = region mask): the board is
# started from the TRUE arm angles (oracle_arm_init.py) with the seed test's schedule; if the fixed point keeps the
# arm (readout score stays low) the failure is search/initialisation, if it drifts to ~1 the per-image likelihood
# cannot hold the arm (SNR floor) and only pooled evidence can.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; REG=${REG:?}; RT=$(basename $REG .mrc); O=$W/oracle_$MODE/${RT}${TAG:-}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') oracle $MODE $RT] $*"; }
GT=$W/sim$MODE/gt_particles.txt; NP=$(grep -vc '^#' $GT); export FB=$A/rnasep_klem/f692e/klem2_bin SIMPLE_QSYS=local
N1=${N1:-7}; N2=${N2:-3}
PS=$A/rnasep_klem/synth_$MODE/src.simple
if [ "${PROJ:-}" = truth ]; then   # PROJ=truth: a project carrying the simulation's TRUE poses (halves by row parity in the board)
  PJ=$W/proj${MODE}_truth; if [ ! -s $PJ/proj.simple ]; then mkdir -p $PJ && cd $PJ && simple_exec prg=new_project projname=proj > np.log 2>&1; cd $PJ/proj
    sed -E 's/ state=[^ ]*//' $W/sim$MODE/all.txt | sed -E 's/$/ state=1/' > $W/sim$MODE/all_import_truth.txt
    simple_exec prg=import_particles projfile=proj.simple stk=$W/sim$MODE/all.mrcs oritab=$W/sim$MODE/all_import_truth.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "truth import FAILED: $(grep -a -m2 ERROR import.log | cut -c1-160)"; exit 1; }
    mv proj.simple ../proj.simple; cd $O; fi
  PS=$PJ/proj.simple; log "truth-pose project $PS"; fi
echo "== $MODE / $RT: truth-initialised lattice"; $PY $W/oracle_arm_init.py $GT $W/truth_$MODE.txt $REG $W/maps160/synth_$MODE $O --n1 $N1 --n2 $N2 ${ARM:+--arm $ARM} 2>&1 | grep -v Warning
echo "== truth readout of the oracle start labels (what a board that keeps the truth would score)"; $PY $W/label_truth_readout.py $O/init_labels.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -9
$PY -c "import numpy as np; np.savetxt('$O/subset.txt', np.arange(1, $NP + 1), fmt='%d')"
log "region board from the oracle start (${N1}x${N2} lattice + 1 free, focus $RT)"; t0=$(date +%s)
[ -s $O/board/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=$N1 N2=$N2 NFREE=${NFREE:-1} CTFABS=${CTFABS:-1} PSCALE=${PSCALE:-1} FOCUSMASK=${FOCUSMASK-$REG} INIT=file INITFILE=$O/init_positions.txt bash $W/fboard.sh $O/board $PS $O/subset.txt "${LPS:-20 14 12 10 10 8 8 8}" "${SIGS:-1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7}" > $O/board.out 2>&1 || { log "board FAILED: $(tail -2 $O/board.out | cut -c1-160)"; exit 1; }
log "board done in $(( $(date +%s) - t0 )) s"; grep -a "argmax changed" $O/board/board.log | sed 's/.*argmax changed/argmax changed/' | cut -c1-60 | tr '\n' ';'; echo
$PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); lab = np.zeros($NP, int); lab[pind - 1] = am; np.savetxt('$O/cells.txt', lab, fmt='%d')
PYEOF
$PY $W/oracle_agree.py $O 2>&1 | grep -v Warning   # start agreement in (g1, g2): init_labels use another cell order
echo "== truth readout of the board's argmax cells after the oracle start"; $PY $W/label_truth_readout.py $O/cells.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -9
log DONE

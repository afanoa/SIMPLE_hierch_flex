#!/bin/bash
#SBATCH -J synth_seed
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# the unified-route test on a synthetic set (env MODE, REG = region mask file): (1) 2D-route seed of the region
# (per-bin principal contrast in the projected mask, sign sync, 3D mode, particle coordinate) read out against the truth
# directly; (2) the Fortran region board on all particles with FOCUSMASK = the region, INIT=file from the seed vs INIT=random,
# same schedule and seed; per-arm / per-site readout of the argmax cells of both.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; REG=${REG:?}; RT=$(basename $REG .mrc); D=$W/binsroute_$MODE; O=$W/seed_$MODE/${RT}${TAG:-}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') seed $MODE $RT] $*"; }
GT=$W/sim$MODE/gt_particles.txt; NP=$(grep -vc '^#' $GT); CONS=$W/maps160/synth_$MODE/consensus_synth.mrc; export FB=$A/rnasep_klem/f692e/klem2_bin SIMPLE_QSYS=local
echo "== $MODE / $RT: 2D-route seed"; $PY $W/region_contrast_seed.py $D/stitch $D/p $D/bins.npz $D/p/cavgs_bins.mrc $REG $CONS $O --npart $NP --n1 ${N1:-7} --n2 ${N2:-1} 2>&1 | grep -v Warning
echo "== truth readout of the seed coordinate (7 bins along the synchronised contrast)"; $PY $W/label_truth_readout.py $O/particle_bin.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -9
$PY -c "import numpy as np; np.savetxt('$O/subset.txt', np.arange(1, $NP + 1), fmt='%d')"
for I in file random; do
  log "region board, init $I (${N1:-7}x${N2:-1} lattice + 1 free, focus $RT)"; [ -s $O/board_$I/klem_latent.txt ] || SERVER_LIST= NPARTS= NTHR=16 N1=${N1:-7} N2=${N2:-1} NFREE=1 CTFABS=1 PSCALE=1 FOCUSMASK=$REG INIT=$I INITFILE=$O/init_positions.txt bash $W/fboard.sh $O/board_$I $A/rnasep_klem/synth_$MODE/src.simple $O/subset.txt "20 14 12 10 10 8 8 8" "1.5 1.2 1.0 1.0 0.8 0.8 0.7 0.7" > $O/board_$I.out 2>&1 || { log "board FAILED: $(tail -2 $O/board_$I.out | cut -c1-160)"; continue; }
  $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board_$I/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); lab = np.zeros($NP, int); lab[pind - 1] = am; np.savetxt('$O/cells_$I.txt', lab, fmt='%d')
PYEOF
  echo "== truth readout of the board's argmax cells, init $I"; $PY $W/label_truth_readout.py $O/cells_$I.txt $GT $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -9
done
log DONE

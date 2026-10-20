#!/bin/bash
#SBATCH -J klem_bag_avgs
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Pooling at the AVERAGE level: the sub-class E-step of the Python board on the expanded averages of the light arms set
# (members inherit the sub-class responsibilities: pool / unmix), arm-focused 3x3 lattice; readout = arm angle std ratio.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE/build; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/synth_arms; O=$A/rnasep_klem/bag_arms; E=$W/route2d_arms/expansion; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') bagavgs] $*"; }
[ -s $O/meta_bags.npz ] || { log "run bag_arms_test.sbatch first (meta_bags.npz)"; exit 1; }
dump(){ simple_exec prg=print_project_field projfile=$E/proj.simple oritype=ptcl2D mkdir=no 2> $O/ppf.err | grep "class=" > $O/bags_full.txt; log "$(grep -c . $O/bags_full.txt) ptcl2D rows with in-plane parameters (build $SIMPLE_PATH)"; [ -s $O/bags_full.txt ]; }
dump || { head -3 $O/ppf.err | cut -c1-200; B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH; dump || { head -3 $O/ppf.err | cut -c1-200; log "no ptcl2D dump"; exit 1; }; }
for MODE in pool unmix; do
  log "board, sub-class E-step on the averages, mode $MODE"; $PY $W/klem.py $W/simarms/all.mrcs none none $O/board_avg$MODE --cons $W/maps160/synth_arms/consensus_synth.mrc --states 1,2,3,4 --axis=0,0,0,0 --init random --meta $O/meta_bags.npz --bags $O/bags_full.txt --bag_parents $O/parents.txt --avgs $E/cls_expansion_cavgs.mrc --avg_mode $MODE --n1 3 --n2 3 --nfree 0 --focus_mask $R/arm1_travel.mrc --lps "12 10 10 8 8 8" --sigmas "1.0 0.8 0.8 0.6 0.5 0.5" --nproc 16 --mskdiam 230 > $O/board_avg$MODE.log 2>&1 || { log "board FAILED: $(grep -a -m2 -i "error" $O/board_avg$MODE.log | cut -c1-160)"; tail -3 $O/board_avg$MODE.log; continue; }
  grep -a "sub-class E-step\|in-plane" $O/board_avg$MODE.log | tail -3 | cut -c1-200
  $PY - <<PYEOF
import numpy as np
L = np.loadtxt('$O/board_avg$MODE/klem_latent.txt', comments='#'); pind = L[:, 0].astype(int); am = L[:, 3].astype(int); N = int(open('$W/simarms/gt_particles.txt').read().count('\n')) + 2
lab = np.zeros(N, int); lab[pind - 1] = am; np.savetxt('$O/cells_avg$MODE.txt', lab, fmt='%d')
PYEOF
  echo "== arm readout, averages mode $MODE (per cell std ratio of each arm; 1 = no information)"; $PY $W/label_truth_readout.py $O/cells_avg$MODE.txt $W/simarms/gt_particles.txt $W/truth_arms.txt arms 2>&1 | grep -v Warning | tail -11
done
log DONE

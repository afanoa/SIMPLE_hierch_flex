#!/bin/bash
#SBATCH -J klem_core_dock
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# the core-aligned frame drifted from the docked frame (10 A corr with the toy consensus 0.930 vs 0.983 for the family
# poses): dock every core-aligned leaf state onto the toy consensus (rigid, dock_volpair) and redo the deposited-conformer
# readouts on the docked set, next to the undocked one
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692e; C=$R/consolidated_12_core_region; O=$C/docked; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') coredock] $*"; }
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
a = read_mrc('$W/maps160/toy_consensus.mrc').astype(np.float32); b = read_mrc('$R/core_frame/fam1_corealigned_docked.mrc').astype(np.float32); c = read_mrc('$R/core_frame/fam1_corealigned.mrc').astype(np.float32)
n = a.shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2; lp = lambda v: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (1.464 / 10.) ** 2)))
def cc(u, v): u, v = lp(u), lp(v); u -= u.mean(); v -= v.mean(); return float((u * v).sum() / np.sqrt((u * u).sum() * (v * v).sum()))
print('family-1 core-aligned map vs toy consensus at 10 A: undocked %.3f, docked %.3f' % (cc(a, c), cc(a, b)))
PYEOF
for m in $C/fam0/regionall/level1/states_dir/state_*.mrc; do s=$(basename $m); [ -s $O/$s ] && continue
  simple_exec prg=dock_volpair vol1=$W/maps160/toy_consensus.mrc vol2=$m smpd=1.464 lpstart=20 lpstop=8 mskdiam=220 nthr=16 mkdir=no > $O/dock_$s.log 2>&1 && d=$(grep -a -o "DOCKED VOLUME WRITTEN: .*" $O/dock_$s.log | awk '{print $4}') && [ -n "$d" ] && mv -f $d $O/$s || log "dock failed for $s"
done
log "docked $(ls $O/state_*.mrc | wc -l) state maps"
for tag in undocked docked; do D=$C/fam0/regionall/level1/states_dir; [ $tag = docked ] && D=$O; S=$O/readout55_$tag; rm -rf "${S:?}"; mkdir -p $S/states; for m in $D/state_*.mrc; do ln -sf $m $S/states/$(basename $m); done; ls $S/states > $S/map_order.txt
  echo "== $tag"; $PY $W/deposited_match.py $S/states $W/maps160/docked_all $W/maps160/toy_consensus.mrc --save $S/matches.txt 2>&1 | grep -v Warning | grep identified; $PY $W/recover55.py $S/matches.txt --names $S/map_order.txt 2>&1 | grep -v Warning | head -1
  $PY $W/deposited_match.py $S/states $W/maps160/docked $W/maps160/toy_consensus.mrc --save $S/m10.txt 2>&1 | grep -v Warning | grep identified; done
log DONE

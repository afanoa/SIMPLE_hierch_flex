#!/bin/bash
#SBATCH -J toyah_sim
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# RNase P apo/holo x open/closed toy ("ah"): the 10 deposited apo conformers (arm swing + 2 closed) and 8 deposited holo
# conformers (RNA + rnpA; closed/open decided by the P12/P10.1 contact of the docked map) simulated together with real
# poses + CTF (as toy10c: calibrated bfac/snr, blob + broken junk), holo share HFRAC of the conformer particles, merged,
# imported, refined from perturbed poses (realistic pose errors). Truth: simah/all.txt, gt_particles.txt (conformer
# >= 101 = holo), closed_ids.txt. env: N (200000), JF (0.10), NP (10), HFRAC (0.5), ANGERR (10), MAXITS (4)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; N=${N:-200000}; JF=${JF:-0.10}; NP=${NP:-10}; HFRAC=${HFRAC:-0.5}
REAL="0.404 0.397 0.157 0.082 0.051"; log(){ echo "[$(date '+%F %T') ah] $*"; }
S=$W/simah; PJ=projah; CSV=$W/rnasep_conformers_ah.csv; mkdir -p $S/calib $W/lat_ah; cd $W
n=0; while [ $n -lt 240 ] && [ ! -s $W/maps160/holo/rnpa_mask.mrc ]; do sleep 60; n=$((n+1)); done
[ -s $W/maps160/holo/rnpa_mask.mrc ] || { log "holo maps not ready"; exit 1; }
# ---- 0 CSV (apo 10 + holo 8, holo share HFRAC) and the closed ids from the docked maps' contact
$PY - <<PYEOF
import csv, glob, sys, numpy as np
sys.path.insert(0, '$W'); from model_class import read_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
W = '$W'; apo = list(csv.DictReader(open(W + '/rnasep_conformers10.csv'))); holo = list(csv.DictReader(open(W + '/rnasep_holo8.csv')))
sa = sum(float(r['images']) for r in apo); sh = sum(float(r['images']) for r in holo); f = $HFRAC / (1 - $HFRAC) * sa / sh
for r in holo: r['images'] = str(int(round(float(r['images']) * f)))
with open('$CSV', 'w', newline='') as fo:
    wr = csv.DictWriter(fo, fieldnames=list(apo[0].keys())); wr.writeheader()
    for r in apo + holo: wr.writerow(r)
g = lambda v: gaussian_filter(v, 2.0); C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in (18, 27, 28, 48, 49)}; cons = read_mrc(W + '/maps160/toy_consensus.mrc')
Lc = g(C[18]); Lo = g(np.mean([C[i] for i in (27, 28, 48, 49)], 0)); diff = Lc - Lo; dens = Lc > Lc.mean() + 1.5 * Lc.std()
core = (diff > diff.mean() + 2.5 * diff.std()) & dens; m = binary_dilation(core, iterations=4); mk = gaussian_filter(m.astype(float), 2.0); mk = np.clip(mk / max(mk.max(), 1e-30), 0, 1)
closed = g(C[18] - cons); openm = g(np.mean([C[i] for i in (27, 28, 48, 49)], 0) - cons)
def fc(x, y): return float(np.sum(mk * x * y) / max(np.sqrt(np.sum(mk * x * x) * np.sum(mk * y * y)), 1e-30))
ids = [18, 51]
for r in holo:
    fs = glob.glob(W + '/maps160/holo_docked/conf%d_*.mrc' % int(r['conformer'])); d = g(read_mrc(fs[0]) - cons); s = fc(d, closed) - fc(d, openm)
    print('holo %s: focal score %+.2f -> %s' % (r['conformer'], s, 'CLOSED' if s > 0.1 else 'open'))
    if s > 0.1: ids.append(int(r['conformer']))
open('$S/closed_ids.txt', 'w').write(','.join(str(i) for i in ids) + '\n'); print('closed ids: %s; CSV %s (holo share %.2f)' % (ids, '$CSV', $HFRAC))
PYEOF
CLOSED=$(cat $S/closed_ids.txt); cd $S
# ---- 1 oritabs
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py $A/rnasep_coupled/basis_span_all/meta_all.npz $CSV $S $N $JF 7 || exit 1
# ---- 2 calibration of bfac / snr against the real band correlations (as toy10c)
head -1500 $S/oris_in_conf01.txt > $S/calib/oris_cal.txt; best=""; bd=1e9
for BF in 0 100 200; do for SNR in 0.0003 0.0006 0.0012 0.0025 0.005 0.01; do
  t=$S/calib/b${BF}_s${SNR}; [ -s $t.txt ] || simple_exec prg=simulate_particles pgrp=c1 vol1=$W/maps160/toy_consensus.mrc smpd=1.464 nptcls=1500 snr=$SNR bfac=$BF ctf=yes kv=300 cs=2.7 fraca=0.1 \
      oritab=$S/calib/oris_cal.txt outstk=$t.mrcs outfile=$t.txt mskdiam=200 nthr=8 mkdir=no > $t.log 2>&1 || { log "calib sim FAILED $BF $SNR"; continue; }
  line=$($PY $W/calib_snr.py $t.mrcs $t.txt $W/maps160/toy_consensus.mrc | grep BANDCORR | cut -d" " -f1-6)
  d=$(echo "$line $REAL" | awk '{d=0; for(i=2;i<=6;i++){ a=$i; b=$(i+5); if(a<0.005)a=0.005; d+=(log(a)-log(b))^2 } printf "%.6f", d}')
  if awk -v a="$d" -v b="$bd" 'BEGIN{exit !(a < b)}'; then bd=$d; best="$BF $SNR"; fi
done; done
[ -n "$best" ] || best="0 0.0012"; set -- $best; BFAC=$1; SNR=$2; log "calibration: bfac $BFAC snr $SNR (log-distance $bd)"
# ---- 3 junk volumes
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
V = read_mrc('$W/maps160/toy_consensus.mrc'); B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx - B/2, yy - B/2), zz - B/2)
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= B * 1.464 / 35.))))); write_mrc('$S/vol_blob.mrc', blob.astype(np.float32), 1.464)
half = 0.5 * (1 + np.tanh((xx - B/2 + 6) / 4.)); write_mrc('$S/vol_broken.mrc', (V * half).astype(np.float32), 1.464); print('junk volumes written')
PYEOF
# ---- 4 simulation per class (4 at a time)
sim(){ name=$1; vol=$2; o=$S/oris_in_$name.txt; [ -s $o ] || return 0; n=$(wc -l < $o); [ -s $S/stk_$name.mrcs ] && return 0; [ $n -gt 0 ] || return 0
  if [ $((n % 2)) = 1 ]; then head -n $((n-1)) $o > ${o%.txt}_even.txt; o=${o%.txt}_even.txt; n=$((n-1)); fi
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$n snr=$SNR bfac=$BFAC ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=200 nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "sim FAILED $name"; }
export -f sim; export S SNR BFAC
{ for m in $W/maps160/docked/conf??_*_m160.mrc $W/maps160/holo_docked/conf1??_*_m160.mrc; do b=$(basename $m .mrc); echo "${b%%_*} $m"; done; echo "blob $S/vol_blob.mrc"; echo "broken $S/vol_broken.mrc"; } | xargs -P 4 -L 1 bash -c 'sim $0 $1'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
# ---- 5 merge, project, refinement from perturbed poses
[ -s $S/all.mrcs ] || $PY $W/merge55.py $S $CSV $CLOSED || exit 1
cd $W; [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > np_$PJ.log 2>&1; cd $W/$PJ
simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_perturbed.txt angerr=${ANGERR:-10} sherr=3 mkdir=no > oriops.log 2>&1 || { log "oriops FAILED"; tail -3 oriops.log; exit 1; }
sed -E 's/ state=[^ ]*//' $S/all_perturbed.txt | sed -E 's/$/ state=1/' > $S/all_import.txt
simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; tail -3 import.log; exit 1; }
simple_exec prg=update_project projfile=$PJ.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
log "project: $(simple_exec prg=print_project_field projfile=$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null | wc -l) particles; refine3D from the toy consensus, perturbed start (angerr ${ANGERR:-10})"
simple_exec prg=refine3D projfile=$PJ.simple vol1=$W/maps160/toy_consensus.mrc nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=200 lp=8 maxits=${MAXITS:-4} nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D FAILED"; grep -a -m3 "ERROR\|THROW" refine.log; exit 1; }
RD=$(ls -td $W/$PJ/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $W/lat_ah/ptcl3D_r0.txt
$PY $W/pose_error.py $W/lat_ah/ptcl3D_r0.txt $S/all.txt --tag "refined vs true" 2>&1 | tail -2
log "DONE toy ah: $S (truth all.txt, gt_particles.txt, closed_ids.txt), project $RD/$PJ.simple, poses lat_ah/ptcl3D_r0.txt"

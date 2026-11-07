#!/bin/bash
#SBATCH -J swingbreath_sim
#SBATCH -p csbdevel
#SBATCH -c 64
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# synthetic benchmark for the organization and the merges, two motions: composition (deposited apo conformers without the closed ones
# and without the breathing base 27, 8 holo), the closed swing of maps160/swing30 (make_swing_maps.py) and a BREATHING family (conformer 27
# with its periphery contracting smoothly, make_breath_maps.py, conformer 300+k = breathing step k); swing
# population densest mid-path (ends ENDW of the middle); real poses + CTF, toy10c calibration (bfac 0, snr 0.0012), blob + broken
# junk, poses = truth perturbed by FASTERR deg. Truth: sim$TAG/gt_particles.txt (conformer 200+k = swing step k, >= 101 holo).
# env TAG (swingbreath), N (400000), JF (0.05 per junk class), SFRAC (0.25), BFRAC (0.20), HFRAC (0.20), BSTEPS (11), AMAX (0.10), ENDW (0.35)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; TAG=${TAG:-swingbreath}; S=$W/sim$TAG; PJ=proj$TAG; CSV=$W/rnasep_conformers_$TAG.csv
N=${N:-400000}; JF=${JF:-0.05}; log(){ echo "[$(date '+%F %T') $TAG] $*"; }; mkdir -p $S/maps $W/lat_$TAG; cd $W
# ---- 0 swing maps from the deposited closed state, and the map set (apo without the closed conformers, holo, swing)
SW=$W/maps160/swing30; BR=$W/maps160/breath${BSTEPS:-11}; [ -s $BR/truth_breath.txt ] || $PY $W/make_breath_maps.py $W/maps160/ah_norm/conf27_emd74230_m160.mrc $BR --amax ${AMAX:-0.10} --nsteps ${BSTEPS:-11} 2>&1 | grep -v Warning
[ -s $SW/truth_swing.txt ] || $PY $W/make_swing_maps.py $W/maps160/ah_norm/conf18_emd74209_m160.mrc $W/maps160/variants/domain_mask.mrc $SW \
  --open $(ls $W/maps160/ah_norm/conf27_*_m160.mrc $W/maps160/ah_norm/conf28_*_m160.mrc $W/maps160/ah_norm/conf48_*_m160.mrc $W/maps160/ah_norm/conf49_*_m160.mrc | tr '\n' ',' | sed 's/,$//') \
  --travel ${TRAVEL:-30} --nsteps ${NSTEPS:-31} 2>&1 | grep -v Warning
for f in $W/maps160/ah_norm/conf*_m160.mrc; do c=$(basename $f | sed 's/^conf//; s/_.*//'); [ "$c" = 18 ] || [ "$c" = 51 ] || [ "$c" = 27 ] || ln -sf $f $S/maps/; done; ln -sf $SW/conf2*_swing_m160.mrc $BR/conf3*_breath_m160.mrc $S/maps/
$PY - <<PYEOF
import csv
W = '$W'; apo = [r for r in csv.DictReader(open(W + '/rnasep_conformers10.csv')) if int(r['conformer']) not in (18, 51, 27)]; holo = list(csv.DictReader(open(W + '/rnasep_holo8.csv')))
ns = ${NSTEPS:-31}; mid = (ns - 1) / 2.; sw = [${ENDW:-0.35} + (1 - ${ENDW:-0.35}) * (1 - abs(k - mid) / mid) for k in range(ns)]
fs, fh, fb = ${SFRAC:-0.25}, ${HFRAC:-0.20}, ${BFRAC:-0.20}; fa = 1 - fs - fh - fb; tot = 1e6; nb = ${BSTEPS:-11}; bm = (nb - 1) / 2.; bw = [${ENDW:-0.35} + (1 - ${ENDW:-0.35}) * (1 - abs(k - bm) / bm) for k in range(nb)]
for rows, frac, wts in ((apo, fa, [float(r['images']) for r in apo]), (holo, fh, [float(r['images']) for r in holo])):
    s = sum(wts)
    for r, w_ in zip(rows, wts): r['images'] = str(int(round(tot * frac * w_ / s)))
swr = [dict(apo[0], emdb_id='swing%02d' % k, conformer=str(200 + k), images=str(int(round(tot * fs * sw[k] / sum(sw)))), file='swing') for k in range(ns)]
brr = [dict(apo[0], emdb_id='breath%02d' % k, conformer=str(300 + k), images=str(int(round(tot * fb * bw[k] / sum(bw)))), file='breath') for k in range(nb)]
with open('$CSV', 'w', newline='') as fo:
    wr = csv.DictWriter(fo, fieldnames=list(apo[0].keys())); wr.writeheader(); [wr.writerow(r) for r in apo + holo + swr + brr]
print('CSV $CSV: %d apo, %d holo, %d swing steps (ends %s, middle %s per 1e6), %d breathing steps (ends %s, middle %s)' % (len(apo), len(holo), ns, swr[0]['images'], swr[ns // 2]['images'], nb, brr[0]['images'], brr[nb // 2]['images']))
PYEOF
CLOSED=$(seq -s, 200 $((200 + ${NSTEPS:-31} - 1)))
# ---- 1 oritabs (real poses + CTF), junk volumes
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py $A/rnasep_coupled/basis_span_all/meta_all.npz $CSV $S $N $JF 7 || exit 1
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
V = read_mrc('$W/maps160/toy_consensus.mrc'); B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx - B/2, yy - B/2), zz - B/2)
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= B * 1.464 / 35.))))); write_mrc('$S/vol_blob.mrc', blob.astype(np.float32), 1.464)
half = 0.5 * (1 + np.tanh((xx - B/2 + 6) / 4.)); write_mrc('$S/vol_broken.mrc', (V * half).astype(np.float32), 1.464); print('junk volumes written')
PYEOF
# ---- 2 simulation per class (8 at a time x 8 threads)
SNR=${SNR:-0.0012}; BFAC=${BFAC:-0}
sim(){ name=$1; vol=$2; o=$S/oris_in_$name.txt; [ -s $o ] || return 0; n=$(wc -l < $o); [ -s $S/stk_$name.mrcs ] && return 0; [ $n -gt 0 ] || return 0
  if [ $((n % 2)) = 1 ]; then head -n $((n-1)) $o > ${o%.txt}_even.txt; o=${o%.txt}_even.txt; n=$((n-1)); fi
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$n snr=$SNR bfac=$BFAC ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=200 nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "sim FAILED $name"; }
export -f sim; export S SNR BFAC
{ for m in $S/maps/conf*_m160.mrc; do b=$(basename $m .mrc); echo "${b%%_*} $m"; done; echo "blob $S/vol_blob.mrc"; echo "broken $S/vol_broken.mrc"; } | xargs -P 8 -L 1 bash -c 'sim $0 $1'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
# ---- 3 merge, project with near-true poses (truth perturbed by FASTERR deg)
[ -s $S/all.mrcs ] || $PY $W/merge55.py $S $CSV $CLOSED || exit 1
cd $W; [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > np_$PJ.log 2>&1; cd $W/$PJ
simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_fast.txt angerr=${FASTERR:-3} sherr=1 mkdir=no > oriops_fast.log 2>&1 || { log "oriops FAILED"; exit 1; }
sed -E 's/ state=[^ ]*//' $S/all_fast.txt | sed -E 's/$/ state=1/' > $S/all_fast_import.txt
simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_fast_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; tail -3 import.log; exit 1; }
simple_exec prg=print_project_field projfile=$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $W/lat_$TAG/ptcl3D_r0.txt
log "DONE: $(grep -c . $W/lat_$TAG/ptcl3D_r0.txt) particles, project $W/$PJ/$PJ.simple, poses lat_$TAG/ptcl3D_r0.txt, truth $S/gt_particles.txt, swing truth $SW/truth_swing.txt, breathing truth $BR/truth_breath.txt"

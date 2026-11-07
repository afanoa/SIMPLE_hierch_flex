#!/bin/bash
#SBATCH -J toy55_sim
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# RNase P 55-conformer toy: waits for the processed maps, calibrates bfac/snr against the real band profile, simulates
# N particles (real poses + CTF, paper image counts as conformer weights, blob + broken junk), merges, imports, refines
# from perturbed poses (realistic pose errors), reconstructs the consensus at the refined poses.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; N=${N:-200000}; JF=${JF:-0.10}; NP=${NP:-10}
CLOSED=18,19,30,50,51,52,53; REAL="0.404 0.397 0.157 0.082 0.051"
log(){ echo "[$(date '+%F %T')] $*"; }
cd $W; n=0; while [ $n -lt 200 ] && [ ! -s $W/maps160/toy_consensus.mrc ]; do sleep 60; n=$((n+1)); done
[ -s $W/maps160/toy_consensus.mrc ] || { log "no toy consensus"; exit 1; }
S=$W/sim; mkdir -p $S/calib && cd $S
# ---- 1 oritabs
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py $A/rnasep_coupled/basis_span_all/meta_all.npz $W/${CSV:-rnasep_conformers10.csv} $S $N $JF 7 || exit 1
# ---- 2 calibration of bfac / snr against the real band correlations (1500 particles from the toy consensus, true poses)
head -1500 $S/oris_in_conf01.txt > $S/calib/oris_cal.txt
best=""; bd=1e9
for BF in 0 100 200; do for SNR in 0.0003 0.0006 0.0012 0.0025 0.005 0.01; do
  t=$S/calib/b${BF}_s${SNR}; [ -s $t.txt ] || simple_exec prg=simulate_particles pgrp=c1 vol1=$W/maps160/toy_consensus.mrc smpd=1.464 nptcls=1500 snr=$SNR bfac=$BF ctf=yes kv=300 cs=2.7 fraca=0.1 \
      oritab=$S/calib/oris_cal.txt outstk=$t.mrcs outfile=$t.txt mskdiam=200 nthr=8 mkdir=no > $t.log 2>&1 || { log "calib sim FAILED $BF $SNR"; continue; }
  line=$($PY $W/calib_snr.py $t.mrcs $t.txt $W/maps160/toy_consensus.mrc | grep BANDCORR | cut -d" " -f1-6); echo "bfac $BF snr $SNR: $line"
  d=$(echo "$line $REAL" | awk '{d=0; for(i=2;i<=6;i++){ a=$i; b=$(i+5); if(a<0.005)a=0.005; d+=(log(a)-log(b))^2 } printf "%.6f", d}'); echo "   log-distance $d"
  if awk -v a="$d" -v b="$bd" 'BEGIN{exit !(a < b)}'; then bd=$d; best="$BF $SNR"; fi
done; done
[ -n "$best" ] || { log "calibration found no best combination; falling back to bfac 0 snr 0.0012"; best="0 0.0012"; }
set -- $best; BFAC=$1; SNR=$2; log "calibration: bfac $BFAC snr $SNR (log-distance $bd); real profile $REAL"
# ---- 3 junk volumes
$PY - <<PYEOF
import numpy as np
def rd(fn):
    h = np.fromfile(fn, np.int32, 256); nx, ny, nz = h[:3]; return np.fromfile(fn, np.float32, offset=1024 + int(h[23]), count=nx*ny*nz).reshape(nz, ny, nx).astype(np.float64)
def wr(fn, v, smpd=1.464):
    v = np.ascontiguousarray(v, np.float32); nz, ny, nx = v.shape; h = np.zeros(256, np.int32); hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz); hf[10:13] = (nx*smpd, ny*smpd, nz*smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean()); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little'); hf[54] = v.std()
    open(fn, 'wb').write(h.tobytes() + v.tobytes())
V = rd('$W/maps160/toy_consensus.mrc'); B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx - B/2, yy - B/2), zz - B/2)
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= B * 1.464 / 35.))))); wr('$S/vol_blob.mrc', blob)
half = 0.5 * (1 + np.tanh((xx - B/2 + 6) / 4.)); wr('$S/vol_broken.mrc', V * half); print('junk volumes written')
PYEOF
# ---- 4 simulation per class (4 at a time)
sim(){ name=$1; vol=$2; o=$S/oris_in_$name.txt; n=$(wc -l < $o); [ -s $S/stk_$name.mrcs ] && return 0; [ $n -gt 0 ] || return 0
  if [ $((n % 2)) = 1 ]; then head -n $((n-1)) $o > ${o%.txt}_even.txt; o=${o%.txt}_even.txt; n=$((n-1)); fi   # simulate_particles needs an even count (build_refspiral)
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$n snr=$SNR bfac=$BFAC ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=200 nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "sim FAILED $name"; }
export -f sim; export S SNR BFAC
{ for m in $W/maps160/docked/conf*_m160.mrc; do b=$(basename $m .mrc); echo "${b%%_*} $m"; done; echo "blob $S/vol_blob.mrc"; echo "broken $S/vol_broken.mrc"; } | xargs -P 4 -L 1 bash -c 'sim $0 $1'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
# ---- 5 merge, project, refinement from perturbed poses
[ -s $S/all.mrcs ] || $PY $W/merge55.py $S $W/${CSV:-rnasep_conformers10.csv} $CLOSED || exit 1
cd $W; rm -rf $W/proj; simple_exec prg=new_project projname=proj > np.log 2>&1; cd $W/proj
simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_perturbed.txt angerr=${ANGERR:-10} sherr=3 mkdir=no > oriops.log 2>&1 || { log "oriops FAILED"; tail -3 oriops.log; exit 1; }
# every particle enters with state=1 (junk included, as in real data; the truth stays in all.txt / gt_particles.txt): refine3D nstates=1 throws on state > 1
sed -E 's/ state=[^ ]*//' $S/all_perturbed.txt | sed -E 's/$/ state=1/' > $S/all_import.txt
simple_exec prg=import_particles projfile=proj.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; tail -3 import.log; exit 1; }
simple_exec prg=update_project projfile=proj.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
log "project: $(simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null | wc -l) particles; refine3D from the toy consensus, perturbed start (angerr ${ANGERR:-10})"
simple_exec prg=refine3D projfile=proj.simple vol1=$W/maps160/toy_consensus.mrc nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=200 lp=8 maxits=${MAXITS:-4} nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D FAILED"; grep -a -m3 "ERROR\|THROW" refine.log; exit 1; }
RD=$(ls -td $W/proj/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $W/proj/ptcl3D_refined.txt
$PY - <<PYEOF
import numpy as np
def rd(fn, keys):
    return np.array([[float(dict(t.split('=',1) for t in l.split() if '=' in t).get(k, 0)) for k in keys] for l in open(fn) if l.strip() and not l.startswith('#')])
def rotmat(ang, axis):
    r = np.radians(ang); c, s = np.cos(r), np.sin(r); R = np.eye(3)
    if axis == 3: R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, -s, s, c
    else:         R[0, 0], R[0, 2], R[2, 0], R[2, 2] = c, -s, s, c
    return R
def e2m(e): return rotmat(e[2], 3) @ rotmat(e[1], 2) @ rotmat(e[0], 3)
r = rd('$W/proj/ptcl3D_refined.txt', ['e1','e2','e3','x','y','state']); t = rd('$S/all.txt', ['e1','e2','e3','x','y','state']); n = min(len(r), len(t))
ang = np.array([np.degrees(np.arccos(np.clip((np.trace(e2m(r[i,:3]) @ e2m(t[i,:3]).T) - 1) / 2, -1, 1))) for i in range(n)]); good = t[:n, 5] > 0
print('refined vs true poses: angle median %.2f deg (p90 %.2f) over conformer particles; junk median %.1f; shift error median %.2f px' % (np.median(ang[good]), np.percentile(ang[good], 90), np.median(ang[~good]), np.median(np.hypot(r[:n,3]-t[:n,3], r[:n,4]-t[:n,4])[good])))
PYEOF
log "DONE toy55: $S (truth all.txt, gt_particles.txt), project $RD/proj.simple"
sbatch $W/run_toy10_battery.sbatch

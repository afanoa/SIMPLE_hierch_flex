#!/bin/bash
#SBATCH -J toybl_sim
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Binding with NO conformational consequence, toy "bl<KDA>": the 10 docked apo conformers and blobbed copies of each
# (a Gaussian blob of KDA kDa at the rnpA site, nothing else changed; conformers 200 + n), 50/50, 100k particles, real
# poses/CTF/SNR, fast poses (true poses perturbed). Truth: conformer >= 201 = bound. env: KDA (13), N (100000), JF (0.10), FASTERR (3)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; KDA=${KDA:-13}; N=${N:-100000}; JF=${JF:-0.10}; BFAC=0; SNR=0.0012
TAG=bl$KDA; S=$W/sim$TAG; PJ=proj$TAG; CSV=$W/rnasep_conformers_$TAG.csv; MAPS=$W/maps160/blob$KDA; L=$W/lat_$TAG; mkdir -p $S $L $MAPS; log(){ echo "[$(date '+%F %T') $TAG] $*"; }; cd $W
# ---- 0 the blobbed maps + CSV (apo rows + blobbed rows), apo maps linked into the same directory for the simulation loop
[ -s $CSV ] || $PY $W/make_blob_maps.py $W/maps160/docked $W/rnasep_conformers10.csv $MAPS $CSV --kda $KDA --site_mask $W/maps160/holo/rnpa_mask.mrc 2>&1 | grep -v Warning
for m in $W/maps160/docked/conf??_*_m160.mrc; do [ -e $MAPS/$(basename $m) ] || ln -s $m $MAPS/$(basename $m); done
echo "18,51,218,251" > $S/closed_ids.txt; CLOSED=$(cat $S/closed_ids.txt)
# ---- 1 oritabs, junk volumes
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py $A/rnasep_coupled/basis_span_all/meta_all.npz $CSV $S $N $JF 7 || exit 1
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
V = read_mrc('$W/maps160/toy_consensus.mrc'); Bx = V.shape[0]; zz, yy, xx = np.mgrid[:Bx, :Bx, :Bx]; r = np.hypot(np.hypot(xx - Bx/2, yy - Bx/2), zz - Bx/2)
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= Bx * 1.464 / 35.))))); write_mrc('$S/vol_blob.mrc', blob.astype(np.float32), 1.464)
half = 0.5 * (1 + np.tanh((xx - Bx/2 + 6) / 4.)); write_mrc('$S/vol_broken.mrc', (V * half).astype(np.float32), 1.464)
PYEOF
# ---- 2 simulation per class (4 at a time)
sim(){ name=$1; vol=$2; o=$S/oris_in_$name.txt; [ -s $o ] || return 0; n=$(wc -l < $o); [ -s $S/stk_$name.mrcs ] && return 0; [ $n -gt 0 ] || return 0
  if [ $((n % 2)) = 1 ]; then head -n $((n-1)) $o > ${o%.txt}_even.txt; o=${o%.txt}_even.txt; n=$((n-1)); fi
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$n snr=$SNR bfac=$BFAC ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=200 nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "sim FAILED $name"; }
export -f sim; export S SNR BFAC
{ for m in $MAPS/conf*_m160.mrc; do b=$(basename $m .mrc); echo "${b%%_*} $m"; done; echo "blob $S/vol_blob.mrc"; echo "broken $S/vol_broken.mrc"; } | xargs -P 4 -L 1 bash -c 'sim $0 $1'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
# ---- 3 merge, fast poses (true poses perturbed), project
[ -s $S/all.mrcs ] || $PY $W/merge55.py $S $CSV $CLOSED || exit 1
RD=$W/$PJ/1_fastposes; mkdir -p $RD && cd $RD
simple_exec prg=oriops oritab=$S/all.txt outfile=$RD/all_fast.txt angerr=${FASTERR:-3} sherr=1 mkdir=no > oriops.log 2>&1 || { log "oriops FAILED"; exit 1; }
sed -E 's/ state=[^ ]*//' $RD/all_fast.txt | sed -E 's/$/ state=1/' > $RD/all_import.txt
rm -rf fast; simple_exec prg=new_project projname=fast > np.log 2>&1; cd fast; simple_exec prg=import_particles projfile=fast.simple stk=$S/all.mrcs oritab=$RD/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; exit 1; }
simple_exec prg=update_project projfile=fast.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=10 nthr=8 mkdir=no > up.log 2>&1
cp -f fast.simple $RD/$PJ.simple; simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $L/ptcl3D_r0.txt
$PY $W/pose_error.py $L/ptcl3D_r0.txt $S/all.txt --tag "fast poses vs true" 2>&1 | tail -1
log "DONE toy $TAG: $S, project $RD/$PJ.simple, poses $L/ptcl3D_r0.txt, site mask $MAPS/blob_site_mask.mrc"

#!/bin/bash
#SBATCH -J sites_sim
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# simsites (factorial site benchmark) end to end in ONE job: truth maps (make_sites_maps.py: rigid conformer-49 core + 3 disjoint sites,
# 8 combination maps) -> independent per-particle site draws with real poses/CTF (make_oris_sites.py) + 20 % junk -> simulation at the
# RNase P ABSOLUTE noise (calibrated snr 0.0012 at bfac 0 for the toy consensus, x every map's projection power; junk too) with the
# volume mask at MSK (232: the closed-arm tip reaches ~110 A) -> merge (merge55.py: gt_particles conformer = combination id) -> project,
# refine3D from 10-deg / 3-px perturbed poses -> pose error -> half stacks (alternate rows) + meta_rows -> the tree (rflex13s) submitted
# as the last step. env N (200000), NP (8), MSK (232), NOTREE (1 = stop before the tree)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K_=$A/rnasep_klem; N=${N:-200000}; NP=${NP:-8}; MSK=${MSK:-232}
S=$W/simsites; PJ=projsites; CSV=$W/simsites.csv; MAPS=$W/maps160/simsites; L=$W/lat_sites; SNR0=0.0012; BFAC=0; mkdir -p $S $L; cd $W; log(){ echo "[$(date '+%F %T') sites] $*"; }
[ -s $MAPS/consensus_sites.mrc ] || $PY $W/make_sites_maps.py $W/maps160 $MAPS $CSV --mskdiam $MSK || exit 1
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris_sites.py $A/rnasep_coupled/basis_span_all/meta_all.npz $MAPS/combos.txt $S $N 0.10 7 || exit 1
$PY - <<PYEOF || exit 1
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
V = read_mrc('$MAPS/consensus_sites.mrc'); B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx - B/2, yy - B/2), zz - B/2)
soft = lambda diam: 0.5 * (1 + np.cos(np.pi * np.clip((r - (diam / 2 / 1.464 - 4)) / 4., 0, 1))); pv = lambda v: np.mean([v.sum(ax).var() for ax in range(3)])
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= B * 1.464 / 35.))))); write_mrc('$S/vol_blob.mrc', blob.astype(np.float32), 1.464)
half = 0.5 * (1 + np.tanh((xx - B/2 + 6) / 4.)); write_mrc('$S/vol_broken.mrc', (V * half).astype(np.float32), 1.464)
ref = pv(read_mrc('$W/maps160/toy_consensus.mrc') * soft(200.))
with open('$S/snr_junk.txt', 'w') as f:
    for nm, v in (('blob', blob), ('broken', V * half)): f.write('%s %.4f\n' % (nm, pv(v * soft($MSK)) / ref))
print('junk volumes written; projection power vs the toy consensus: ' + open('$S/snr_junk.txt').read().replace('\n', ' '))
PYEOF
sim(){ name=$1; vol=$2; sc=$3; o=$S/oris_in_$name.txt; [ -s $o ] || return 0; n=$(wc -l < $o); [ -s $S/stk_$name.mrcs ] && return 0; [ $n -gt 0 ] || return 0
  if [ $((n % 2)) = 1 ]; then head -n $((n-1)) $o > ${o%.txt}_even.txt; o=${o%.txt}_even.txt; n=$((n-1)); fi
  snr=$(awk -v a=$SNR0 -v b=$sc 'BEGIN {printf "%.7f", a * b}')
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$n snr=$snr bfac=$BFAC ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=$MSK nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "sim FAILED $name"; }
export -f sim; export S SNR0 BFAC MSK
log "simulating at snr $SNR0 x projection power (absolute RNase P noise), bfac $BFAC, volume mask $MSK A"
{ while read nm sc; do echo "$nm $MAPS/${nm}_sites.mrc $sc"; done < $MAPS/snr_scale.txt; while read nm sc; do echo "$nm $S/vol_$nm.mrc $sc"; done < $S/snr_junk.txt; } | xargs -P 3 -L 1 bash -c 'sim $0 $1 $2'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
[ -s $S/all.mrcs ] || $PY $W/merge55.py $S $CSV 5,6,7,8 || exit 1
if [ ! -s $L/ptcl3D_r0.txt ]; then
  cd $W; [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > np_$PJ.log 2>&1; cd $W/$PJ
  simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_perturbed.txt angerr=10 sherr=3 mkdir=no > oriops.log 2>&1 || { log "oriops FAILED"; exit 1; }
  sed -E 's/ state=[^ ]*//' $S/all_perturbed.txt | sed -E 's/$/ state=1/' > $S/all_import.txt
  simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; exit 1; }
  simple_exec prg=update_project projfile=$PJ.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
  log "refine3D (prob, lp 8, mskdiam $MSK, 4 its) from 10-deg perturbed poses, $NP parts"
  simple_exec prg=refine3D projfile=$PJ.simple vol1=$MAPS/consensus_sites.mrc nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=$MSK lp=8 maxits=4 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D FAILED"; grep -a -m3 "ERROR\|THROW" refine.log; exit 1; }
  RD=$(ls -td $W/$PJ/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $L/ptcl3D_r0.txt
fi
echo "== pose error of the refined poses vs the truth"; $PY $W/pose_error.py $L/ptcl3D_r0.txt $S/all.txt --tag "refined vs true" 2>&1 | tail -2
[ -s ${S}_B/all.mrcs ] || $PY $W/split_halves.py $S ${S}_A ${S}_B || exit 1
NR=$(grep -c . $S/all.txt); [ -s ${S}_B/meta_rows.txt ] || { seq 1 2 $NR > ${S}_A/meta_rows.txt; seq 2 2 $NR > ${S}_B/meta_rows.txt; }
SRC=$K_/sites_src; mkdir -p $SRC; RD=$(ls -td $W/$PJ/*_refine3D | head -1)
[ -s $SRC/src.simple ] || { cp -f $RD/$PJ.simple $SRC/src.simple && cd $SRC && simple_exec prg=update_project projfile=$SRC/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=8 nthr=8 mkdir=no > $SRC/up.log 2>&1; }
log "halves ${S}_A / ${S}_B ($(grep -c . ${S}_A/meta_rows.txt) / $(grep -c . ${S}_B/meta_rows.txt) rows), tree source $SRC/src.simple"
[ "${NOTREE:-0}" = 1 ] && { log "NOTREE: done without the tree"; exit 0; }
cd $W; RT=simsites TAG=sites SRC=$SRC/src.simple H=$S CP=$W/benchcache_sites CTFABS=0 MSKDIAM=$MSK KFLOOR=5000 MINP=15000 MAXD=4 BORROW=none PAR=${TPAR:-1} NP=8 \
  sbatch --export=ALL --mem=128G $W/rflex13s.sbatch; log "tree submitted (rflex13s, RT simsites)"
log DONE

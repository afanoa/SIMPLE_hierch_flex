#!/bin/bash
#SBATCH -J qtoy_sim
#SBATCH -p csbdevel
#SBATCH -c 24
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# qtoy end to end in ONE job: N RNase P particles with truth for the states (qtoy.csv: 6 docked deposited conformers along the arm swing, a
# rare closed state, a holo; plus blob / broken junk) AND a per-particle quality spread independent of the state (qtoy_merge.py: B, SNR scale,
# one absolute noise): real poses / CTF (make_oris55.py), noise-free simulate_particles per class, merge, refine3D (prob, lp 8, 4 its) from
# 10-deg / 3-px perturbed truth -> simqtoy/src.simple (the trees' PROJ) and simqtoy/truth_qtoy.txt (TRUTH). env N (60000), NP (8), MSK (232)
# QT (qtoy: sim$QT/, proj$QT, $QT.csv) and SNR0 (qtoy_merge.py's, 0.0012) for toy versions
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel SIMPLE_EMAIL=nil OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; N=${N:-60000}; NP=${NP:-8}; MSK=${MSK:-232}
QT=${QT:-qtoy}; S=$W/sim$QT; PJ=proj$QT; CSV=$W/$QT.csv; M=$W/maps160; export SNR0=${SNR0:-0.0012}; mkdir -p $S; cd $W; log(){ echo "[$(date '+%F %T') $QT] $*"; }
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py $A/rnasep_coupled/basis_span_all/meta_all.npz $CSV $S $N 0.05 7 || exit 1
[ -s $S/vol_broken.mrc ] || $PY - <<PYEOF || exit 1
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
V = read_mrc('$M/toy_consensus.mrc'); B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx - B/2, yy - B/2), zz - B/2)
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= B * 1.464 / 35.)))))
write_mrc('$S/vol_blob.mrc', blob.astype(np.float32), 1.464); half = 0.5 * (1 + np.tanh((xx - B/2 + 6) / 4.)); write_mrc('$S/vol_broken.mrc', (V * half).astype(np.float32), 1.464)
print('junk volumes written (blob: the consensus low-passed to 35 A; broken: its half)')
PYEOF
sim(){ name=$1; vol=$2; o=$S/oris_in_$name.txt; [ -s $S/stk_$name.mrcs ] && return 0
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$(grep -c . $o) snr=10 bfac=0 ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=$MSK nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "simulate_particles FAILED for $name"; }
export -f sim; export S MSK
log "simulating noise-free (snr 10: simulate_particles adds no noise at snr >= 5), bfac 0, volume mask $MSK A"
{ while IFS=, read c img kind s1 ca map; do [ "$c" = conformer ] || echo "conf$(printf %02d $c) $M/$map"; done < $CSV; echo "blob $S/vol_blob.mrc"; echo "broken $S/vol_broken.mrc"; } | xargs -P 3 -L 1 bash -c 'sim $0 $1'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
[ -s $S/all.mrcs ] || $PY $W/qtoy_merge.py $S $CSV $S/truth_qtoy.txt || exit 1
if [ ! -s $S/refined.txt ]; then
  [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > np_$PJ.log 2>&1; cd $W/$PJ
  simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_perturbed.txt angerr=10 sherr=3 mkdir=no > oriops.log 2>&1 || { log "oriops FAILED"; exit 1; }
  sed -E 's/ state=[^ ]*//' $S/all_perturbed.txt | sed -E 's/$/ state=1/' > $S/all_import.txt
  simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; exit 1; }
  simple_exec prg=update_project projfile=$PJ.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
  log "refine3D (prob, lp 8, mskdiam $MSK, 4 its) from 10-deg / 3-px perturbed poses, $NP parts"
  simple_exec prg=refine3D projfile=$PJ.simple vol1=$M/toy_consensus.mrc nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=$MSK lp=8 maxits=4 nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D FAILED"; grep -a -m3 "ERROR\|THROW" refine.log; exit 1; }
  RD=$(ls -td $W/$PJ/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $S/refined.txt
fi
echo "== pose error of the refined poses vs the truth"; $PY $W/pose_error.py $S/refined.txt $S/all.txt --tag "refined vs true" 2>&1 | tail -2
RD=$(ls -td $W/$PJ/*_refine3D | head -1); [ -s $S/src.simple ] || cp $RD/$PJ.simple $S/src.simple
log "DONE: PROJ $S/src.simple MSK $MSK TRUTH $S/truth_qtoy.txt N $(grep -vc '^#' $S/truth_qtoy.txt)"

#!/bin/bash
#SBATCH -J toysynth_sim
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Synthetic benchmark (MODE arms | tethered) through the toy pipeline: maps (make_synth_maps.py) -> real poses/CTF
# oritabs -> bfac/snr from the ah calibration -> simulate per configuration + junk -> merge -> project -> refine3D from
# the synthetic consensus with perturbed poses. Truth: sim$MODE/all.txt, gt_particles.txt, truth_$MODE.txt. env: MODE, N (200000), NP (10)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:-arms}; N=${N:-200000}; JF=${JF:-0.10}; NP=${NP:-10}
S=$W/sim$MODE; PJ=proj$MODE; CSV=$W/synth_$MODE.csv; MAPS=$W/maps160/synth_$MODE; mkdir -p $S $W/lat_$MODE; cd $W; log(){ echo "[$(date '+%F %T') synth $MODE] $*"; }
[ -s $MAPS/consensus_synth.mrc ] || $PY $W/make_synth_maps.py $MODE $MAPS $CSV $W/truth_$MODE.txt $W/maps160/toy_consensus.mrc --nconf 100 ${SYNTH_EXTRA:-} || exit 1
CONS=$MAPS/consensus_synth.mrc; cd $S
[ -s $S/oris_in_blob.txt ] || $PY $W/make_oris55.py $A/rnasep_coupled/basis_span_all/meta_all.npz $CSV $S $N $JF 7 || exit 1
cal=$(grep -a -h "calibration: bfac" $(ls -t $W/rnasep_toyah_sim*.out $W/toyah_sim*.out 2>/dev/null) 2>/dev/null | tail -1); BFAC=$(echo "$cal" | sed -E 's/.*bfac ([0-9.]+) snr ([0-9.]+).*/\1/'); SNR=$(echo "$cal" | sed -E 's/.*bfac ([0-9.]+) snr ([0-9.]+).*/\2/')
[ -n "$BFAC" ] && [ -n "$SNR" ] && [ "$BFAC" != "$cal" ] || { BFAC=0; SNR=0.0012; }; log "bfac $BFAC snr $SNR (from the ah calibration)"
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
V = read_mrc('$CONS'); B = V.shape[0]; zz, yy, xx = np.mgrid[:B, :B, :B]; r = np.hypot(np.hypot(xx - B/2, yy - B/2), zz - B/2)
F = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(V))); blob = np.real(np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(F * (r <= B * 1.464 / 35.))))); write_mrc('$S/vol_blob.mrc', blob.astype(np.float32), 1.464)
half = 0.5 * (1 + np.tanh((xx - B/2 + 6) / 4.)); write_mrc('$S/vol_broken.mrc', (V * half).astype(np.float32), 1.464); print('junk volumes written')
PYEOF
sim(){ name=$1; vol=$2; o=$S/oris_in_$name.txt; [ -s $o ] || return 0; n=$(wc -l < $o); [ -s $S/stk_$name.mrcs ] && return 0; [ $n -gt 0 ] || return 0
  if [ $((n % 2)) = 1 ]; then head -n $((n-1)) $o > ${o%.txt}_even.txt; o=${o%.txt}_even.txt; n=$((n-1)); fi
  simple_exec prg=simulate_particles pgrp=c1 vol1=$vol smpd=1.464 nptcls=$n snr=$SNR bfac=$BFAC ctf=yes kv=300 cs=2.7 fraca=0.1 sherr=2 oritab=$o outstk=$S/stk_$name.mrcs outfile=$S/oris_$name.txt mskdiam=200 nthr=8 mkdir=no > $S/sim_$name.log 2>&1 || echo "sim FAILED $name"; }
export -f sim; export S SNR BFAC
{ for m in $MAPS/conf??_synth.mrc $MAPS/conf???_synth.mrc; do [ -f $m ] || continue; b=$(basename $m .mrc); echo "${b%%_*} $m"; done; echo "blob $S/vol_blob.mrc"; echo "broken $S/vol_broken.mrc"; } | xargs -P 4 -L 1 bash -c 'sim $0 $1'
log "simulated: $(ls $S/stk_*.mrcs | wc -l) stacks"
[ -s $S/all.mrcs ] || $PY $W/merge55.py $S $CSV 0 || exit 1
cd $W; [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > np_$PJ.log 2>&1; cd $W/$PJ
simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_perturbed.txt angerr=${ANGERR:-10} sherr=3 mkdir=no > oriops.log 2>&1 || { log "oriops FAILED"; exit 1; }
sed -E 's/ state=[^ ]*//' $S/all_perturbed.txt | sed -E 's/$/ state=1/' > $S/all_import.txt
simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; exit 1; }
simple_exec prg=update_project projfile=$PJ.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=48000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
simple_exec prg=refine3D projfile=$PJ.simple vol1=$CONS nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=200 lp=8 maxits=${MAXITS:-4} nparts=$NP nthr=8 > refine.log 2>&1 || { log "refine3D FAILED"; grep -a -m3 "ERROR\|THROW" refine.log; exit 1; }
RD=$(ls -td $W/$PJ/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $W/lat_$MODE/ptcl3D_r0.txt
$PY $W/pose_error.py $W/lat_$MODE/ptcl3D_r0.txt $S/all.txt --tag "refined vs true" 2>&1 | tail -2
log "DONE synthetic $MODE: $S, project $RD/$PJ.simple, truth $W/truth_$MODE.txt"

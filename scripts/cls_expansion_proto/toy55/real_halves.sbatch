#!/bin/bash
#SBATCH -J real_halves
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=110G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# real RNase P gold-standard halves for the data-derived tree (MODE rnasep692): row-ordered half stacks + pose tables
# (real_halves.py), per-half SIMPLE projects imported with ctf=flip, run dirs (src.simple, meta_toy.npz, subset_all.txt),
# then the bench cache and the cross-half geometry check (geometry_check.py). One job, resumable step by step.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem; MODE=${MODE:-rnasep692}; export OMP_NUM_THREADS=8
log(){ echo "[$(date '+%F %T') real_halves $MODE] $*"; }
[ -s $W/sim${MODE}_B/meta_rows.txt ] || $PY $W/real_halves.py $A/rnasep_coupled/basis_span_all/meta_all.npz "/mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part*.mrcs" $W/sim$MODE $W/lat_$MODE || { log "halves FAILED"; exit 1; }
for H in A B; do PJ=proj${MODE}_$H; R=$K/synth_${MODE}_$H; mkdir -p $R
  if [ ! -s $R/src.simple ]; then cd $W; [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > np_$PJ.log 2>&1; cd $W/$PJ
    simple_exec prg=import_particles projfile=$PJ.simple stk=$W/sim${MODE}_$H/all.mrcs oritab=$W/lat_${MODE}_$H/ptcl3D_r0.txt ctf=flip smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import $H FAILED: $(grep -a -m2 ERROR import.log | cut -c1-160)"; exit 1; }
    cp -f $W/$PJ/$PJ.simple $R/src.simple; cd $R; simple_exec prg=update_project projfile=$R/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > up.log 2>&1; fi
  [ -s $R/meta_toy.npz ] || $PY $W/make_toy_meta.py $W/lat_${MODE}_$H/ptcl3D_r0.txt $R/meta_toy.npz
  [ -s $R/subset_all.txt ] || $PY -c "import numpy as np; n = len(np.load('$R/meta_toy.npz')['flags']); np.savetxt('$R/subset_all.txt', np.arange(1, n + 1), fmt='%d'); print('$H subset', n)"
  log "half $H: run dir $R ($(grep -c . $R/subset_all.txt) particles, project $(ls -la $R/src.simple | awk '{print $5}') bytes)"; done
$PY $W/geometry_check.py $W/sim$MODE $W/lat_$MODE $W/benchcache_$MODE 40000 2>&1 | grep -v Warning
log DONE

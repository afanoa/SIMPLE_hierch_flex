#!/bin/bash
#SBATCH -J igm_dtree
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The data-derived tree (dtree_v10) on a synthetic set from run_igm_sim, in ONE job: two independent halves (split_halves.py,
# alternate rows) -> each half refined on its own (refine3D prob from 10-deg perturbed poses against the set's consensus,
# in parallel) -> half run dirs (src.simple, meta_toy.npz, subset_all.txt) -> dtree_v10 inline -> per-arm readout of every
# node's children (dtree_arm_readout.py). env MODE (armsigm1), TAG (dtree), MSK (220), NP (10), MINGRP (1000), NMIN (2000), PARNODE (3)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; K=$A/rnasep_klem
MODE=${MODE:-armsigm1}; TAG=${TAG:-dtree}; MSK=${MSK:-220}; NP=${NP:-10}; CONS=$W/maps160/synth_$MODE/consensus_synth.mrc; log(){ echo "[$(date '+%F %T') igm_dtree $MODE] $*"; }
[ -s $W/sim${MODE}_B/all.mrcs ] || $PY $W/split_halves.py $W/sim$MODE $W/sim${MODE}_A $W/sim${MODE}_B || exit 1
half(){ local H=$1 S=$W/sim${MODE}_$1 PJ=proj${MODE}_$1 L=$W/lat_${MODE}_$1 R=$K/synth_${MODE}_$1; mkdir -p $L $R
  if [ ! -s $L/ptcl3D_r0.txt ]; then
    cd $W; [ -e $W/$PJ ] && mv $W/$PJ $W/${PJ}_old_$(date +%s); simple_exec prg=new_project projname=$PJ > $W/np_$PJ.log 2>&1; cd $W/$PJ || return 1
    simple_exec prg=oriops oritab=$S/all.txt outfile=$S/all_perturbed.txt angerr=${ANGERR:-10} sherr=3 mkdir=no > oriops.log 2>&1 || { echo "half $H oriops FAILED"; return 1; }
    sed -E 's/ state=[^ ]*//' $S/all_perturbed.txt | sed -E 's/$/ state=1/' > $S/all_import.txt
    simple_exec prg=import_particles projfile=$PJ.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { echo "half $H import FAILED"; return 1; }
    simple_exec prg=update_project projfile=$PJ.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=43200 nparts=$NP nthr=8 mkdir=no > up.log 2>&1
    simple_exec prg=refine3D projfile=$PJ.simple vol1=$CONS nstates=1 refine=prob objfun=cc pgrp=c1 mskdiam=$MSK lp=8 maxits=${MAXITS:-4} nparts=$NP nthr=8 > refine.log 2>&1 || { echo "half $H refine3D FAILED: $(grep -a -m2 'ERROR\|THROW' refine.log | cut -c1-160)"; return 1; }
    local RD=$(ls -td $W/$PJ/*_refine3D | head -1); simple_exec prg=print_project_field projfile=$RD/$PJ.simple oritype=ptcl3D mkdir=no 2>/dev/null > $L/ptcl3D_r0.txt
  fi
  echo "== half $H: $($PY $W/pose_error.py $L/ptcl3D_r0.txt $S/all.txt --tag "refined vs true" 2>&1 | tail -1)"
  local RD=$(ls -td $W/$PJ/*_refine3D | head -1); $PY $W/make_toy_meta.py $L/ptcl3D_r0.txt $R/meta_toy.npz > /dev/null || return 1; cp -f $RD/$PJ.simple $R/src.simple
  simple_exec prg=update_project projfile=$R/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > $R/up.log 2>&1
  $PY -c "import numpy as np; n = sum(1 for l in open('$L/ptcl3D_r0.txt') if l.strip()); np.savetxt('$R/subset_all.txt', np.arange(1, n + 1), fmt='%d')"; return 0; }
log "halves: split + independent refinements"; half A > $W/igm_half_${MODE}_A.log 2>&1 & pa=$!; half B > $W/igm_half_${MODE}_B.log 2>&1 & pb=$!
wait $pa; ea=$?; wait $pb; eb=$?; cat $W/igm_half_${MODE}_A.log $W/igm_half_${MODE}_B.log | grep -a "==\|FAILED"; [ $ea = 0 ] && [ $eb = 0 ] || { log "half preparation FAILED ($ea $eb)"; exit 1; }
log "data-derived tree (dtree_v10, TAG $TAG)"; MODE=$MODE TAG=$TAG MINGRP=${MINGRP:-1000} NMIN=${NMIN:-2000} PARNODE=${PARNODE:-3} bash $W/dtree_v10.sbatch > $K/dtree_${MODE}_$TAG.out 2>&1; log "dtree exit $? (log $K/dtree_${MODE}_$TAG.out)"
grep -a "^== node\|K\*" $K/dtree_${MODE}_$TAG/TREE.txt 2>/dev/null | head -40
echo "== per-arm readout of the tree's children"; $PY $W/dtree_arm_readout.py $K $MODE $TAG $W/sim$MODE $W/truth_$MODE.txt 2>&1 | grep -v Warning
log DONE

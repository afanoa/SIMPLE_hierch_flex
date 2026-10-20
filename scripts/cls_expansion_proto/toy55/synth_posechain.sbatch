#!/bin/bash
#SBATCH -J synth_posechain
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The full klem tree on a synthetic set at near-true poses: the simulation's TRUE poses perturbed by ANGERR deg (default 3,
# what a large rigid core delivers) -> project + meta -> f692_chain10 inline (RTAG synth_<MODE>_tp<ANGERR>) -> per-arm
# readout, next to the same readout of the tree at the refined poses. Separates what the tree does once poses are good
# from what pose error costs. env MODE (arms3), ANGERR (3), SHERR (1)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; AE=${ANGERR:-3}
TAG=synth_${MODE}_tp$AE; S=$W/sim$MODE; P=$W/projtp_${MODE}_$AE; R=$A/rnasep_klem/$TAG; mkdir -p $P $R; log(){ echo "[$(date '+%F %T') posechain $TAG] $*"; }
if [ ! -s $P/proj.simple ]; then
  cd $P && simple_exec prg=new_project projname=proj > np.log 2>&1; cd $P/proj
  simple_exec prg=oriops oritab=$S/all.txt outfile=$P/all_tp.txt angerr=$AE sherr=${SHERR:-1} mkdir=no > $P/oriops.log 2>&1 || { log "oriops FAILED"; exit 1; }
  sed -E 's/ state=[^ ]*//' $P/all_tp.txt | sed -E 's/$/ state=1/' > $P/all_tp_import.txt
  simple_exec prg=import_particles projfile=proj.simple stk=$S/all.mrcs oritab=$P/all_tp_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > $P/import.log 2>&1 || { log "import FAILED: $(grep -a -m2 ERROR $P/import.log | cut -c1-160)"; exit 1; }
  mv proj.simple $P/proj.simple
fi
simple_exec prg=print_project_field projfile=$P/proj.simple oritype=ptcl3D mkdir=no 2>/dev/null > $P/ptcl3D_tp.txt
echo "== pose error of the project vs the truth"; $PY $W/pose_error.py $P/ptcl3D_tp.txt $S/all.txt --tag "truth + $AE deg" 2>&1 | tail -2
$PY $W/make_toy_meta.py $P/ptcl3D_tp.txt $R/meta_toy.npz || exit 1; cp -f $P/proj.simple $R/src.simple
simple_exec prg=update_project projfile=$R/src.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=24000 walltime=43200 nparts=10 nthr=8 mkdir=no > $R/up.log 2>&1
log "klem tree (f692_chain10) on $R"; RTAG=$TAG SRC=$R/src.simple METAALL=$R/meta_toy.npz bash $W/f692_chain10.sbatch > $R/chain.out 2>&1; log "chain exit $? (log $R/chain.out)"
echo "== per-arm readout, tree at near-true poses ($TAG)"; $PY $W/arm_angle_readout.py $R $S/gt_particles.txt $W/truth_$MODE.txt 2>&1 | grep -v Warning | tail -40
echo "== per-arm readout, tree at the refined poses (synth_$MODE, baseline)"; $PY $W/arm_angle_readout.py $A/rnasep_klem/synth_$MODE $S/gt_particles.txt $W/truth_$MODE.txt 2>&1 | grep -v Warning | tail -40
log DONE

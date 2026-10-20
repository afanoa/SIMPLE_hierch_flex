#!/bin/bash
#SBATCH -J synth_2droute
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# the pooled (2D) route on a synthetic benchmark (env MODE): abinitio2D (100 classes) on the refined project, then
# cls_expansion with kernel sub-class averages on the 5 best classes (run_kernel_cls.sh), then the per-particle sub-class
# labels against the truth (arm angles / site occupancies): does pooling before the decision see what the boards cannot?
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE/build; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8   # abinitio2D is in the current tree, not in the kernel snapshot
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; MODE=${MODE:?}; D=$W/route2d_$MODE; mkdir -p $D; cd $D; log(){ echo "[$(date '+%F %T') 2droute $MODE] $*"; }
# a fresh project from the simulated stack with THIS build (the kernel-snapshot project's record width is unreadable here)
S=$W/sim$MODE; if [ -f $D/p/1_abinitio2D/ABINITIO2D_FINISHED ]; then log "2D classes exist, reusing"; cd $D/p; else simple_exec prg=new_project projname=p > np.log 2>&1; cd $D/p
simple_exec prg=import_particles projfile=p.simple stk=$S/all.mrcs oritab=$S/all_import.txt ctf=yes smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED: $(grep -a -m2 ERROR import.log | cut -c1-160)"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=43200 nparts=${NP:-10} nthr=8 mkdir=no > up.log 2>&1
log "abinitio2D ncls ${NCLS:-100} on $(grep -c . $W/sim$MODE/gt_particles.txt) particles"
simple_exec prg=abinitio2D projfile=p.simple ncls=${NCLS:-100} mskdiam=${MSK:-200} nparts=${NP:-10} nthr=8 > abinitio2D.log 2>&1 || { log "abinitio2D FAILED: $(grep -a -m2 ERROR abinitio2D.log | cut -c1-160)"; exit 1; }
fi
PD=$(ls -td $D/p/*_abinitio2D | head -1); log "2D done: $PD"
rm -rf $D/expansion; export CLASSES=$D/expansion/classes.txt   # run_kernel_cls.sh reads the class list from $CLASSES
log "cls_expansion (kernel sub-class averages, 5 best classes)"; bash $A/clsexp_stage1/run_kernel_cls.sh $D/expansion $PD/p.simple ${MSK:-200} 2>&1 | grep -v "^$" | tail -8
simple_exec prg=print_project_field projfile=$D/expansion/proj.simple oritype=ptcl2D mkdir=no 2>/dev/null | awk '{c=0; for(i=1;i<=NF;i++) if($i ~ /^class=/){split($i,a,"="); c=a[2]}; print int(c)}' > $D/subclass_labels.txt
log "labels: $(sort -n $D/subclass_labels.txt | uniq | wc -l) distinct sub-class ids over $(wc -l < $D/subclass_labels.txt) rows"
$PY $W/label_truth_readout.py $D/subclass_labels.txt $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -25
log "DONE"

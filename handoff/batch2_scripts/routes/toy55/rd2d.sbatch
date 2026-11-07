#!/bin/bash
#SBATCH -J rd2d
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# sub-class labels of a 2D-route expansion project (route2d_<MODE>) and their per-arm / per-site truth readout (env MODE)
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${BLD:-$HOME/SIMPLE_stage1/build}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; MODE=${MODE:?}; D=$W/route2d_$MODE
simple_exec prg=print_project_field projfile=$D/expansion/proj.simple oritype=ptcl2D mkdir=no > $D/ptcl2D_dump.txt 2> $D/ptcl2D_dump.err
echo "dump: $(grep -c . $D/ptcl2D_dump.txt) lines; stderr $(head -c 300 $D/ptcl2D_dump.err)"; head -c 400 $D/ptcl2D_dump.txt; echo
awk '/class=/{c=0; for(i=1;i<=NF;i++) if($i ~ /^class=/){split($i,a,"="); c=a[2]}; print int(c)}' $D/ptcl2D_dump.txt > $D/subclass_labels.txt
echo "labels: $(sort -n $D/subclass_labels.txt | uniq | wc -l) distinct ids over $(wc -l < $D/subclass_labels.txt) rows"
$PY $W/label_truth_readout.py $D/subclass_labels.txt $W/sim$MODE/gt_particles.txt $W/truth_$MODE.txt $MODE 2>&1 | grep -v Warning | tail -25
echo RD2D DONE

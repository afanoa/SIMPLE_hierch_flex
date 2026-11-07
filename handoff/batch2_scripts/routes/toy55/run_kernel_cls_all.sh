#!/bin/bash
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 3:00:00
#SBATCH --mail-type=NONE
# run_kernel_cls_all.sh OUTDIR SRCPROJ MSKDIAM: stage-1 cls_expansion (conditional averages) on ALL selected classes with
# >= MINPOP particles (300), sub-class budget ncls = PERCLS (20) x classes, nparts=5 x nthr=8, from ~/SIMPLE_stage1/build
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64:${LD_LIBRARY_PATH:-}
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*"); do unset $v; done
B=${BUILD:-$HOME/SIMPLE_stage1/build}
export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel SIMPLE_IO_NSTREAMS=32
D=$1; SRC=$2; MSK=$3
QS="qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 nparts=5 nthr=8"
mkdir -p $D && cd $D || exit 4
[ -e proj.simple ] && { echo "REFUSE: $D has a project"; exit 3; }
cp $SRC proj.simple
simple_exec prg=update_project projfile=$D/proj.simple $QS mkdir=no > update_project.log 2>&1 || { echo "update_project FAILED"; exit 5; }
simple_exec prg=print_project_field projfile=$D/proj.simple oritype=cls2D > cls2D_in.txt 2>/dev/null
python3 -c "
import re
rows = []
for i, l in enumerate(open(\"cls2D_in.txt\")):
    f = dict(re.findall(r\"(\w+)=\s*([-0-9.Ee+]+)\", l))
    if \"res\" in f and float(f.get(\"state\", 1)) > 0 and float(f.get(\"pop\", 0)) >= ${MINPOP:-300}:
        rows.append((float(f[\"res\"]), int(float(f.get(\"class\", i + 1))), int(float(f[\"pop\"]))))
rows.sort()
open(\"classes.txt\", \"w\").write(\"\".join(\"%d\n\" % c for _, c, _ in rows)); open(\"ncls.txt\", \"w\").write(\"%d\n\" % (${PERCLS:-20} * len(rows)))
print(\"all\", len(rows), \"eligible classes, sub-class budget\", ${PERCLS:-20} * len(rows))
"
echo "build=$(cd $HOME/SIMPLE_stage1 && /usr/bin/git log --oneline -1 | cut -c1-9) + stage-1 diff; host=$(hostname)"; date "+START %F %T"
simple_exec prg=cls_expansion projfile=$D/proj.simple ncls=$(cat ncls.txt) mskdiam=$MSK class_assignment=$CLASSES subcls_avg=kernel $QS mkdir=no > CLS_EXPANSION_OUTPUT 2>&1
echo "RC=$?"; date "+END %F %T"
grep -a "ERROR" CLS_EXPANSION_OUTPUT SIMPLE_SUBPROC_OUTPUT 2>/dev/null | head -3
grep -ah "SUBCLASS BUDGET\|kept whole\|split threshold" CLS_EXPANSION_OUTPUT | head -3
grep -ah "Cls expansion bands\|Cls expansion conditional averages\|class kept whole: class" SIMPLE_SUBPROC_OUTPUT CLS_EXPANSION_OUTPUT 2>/dev/null
simple_exec prg=print_project_field projfile=$D/proj.simple oritype=out 2>/dev/null | grep -o "imgkind=[A-Za-z0-9_]*"
ls -la cls_expansion_cavgs*.mrc cls_expansion_frcs.bin 2>/dev/null

#!/bin/bash
#SBATCH -J klem_dump2d
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
B=$HOME/SIMPLE_stage1/build; export SIMPLE_PATH=$B PATH=$B/bin:$PATH
for M in arms tethered; do E=/mnt/beegfs/elmlund/afan/routes/toy55/route2d_$M/expansion; echo "== $M"; simple_exec prg=print_project_field projfile=$E/proj.simple oritype=ptcl2D mkdir=no 2>/dev/null > $E/ptcl2D.txt; head -2 $E/ptcl2D.txt | cut -c1-200; echo "rows $(wc -l < $E/ptcl2D.txt)"; awk "{for(i=1;i<=NF;i++){if(\$i ~ /^class=/){split(\$i,a,\"=\"); c=int(a[2])}; if(\$i ~ /^cluster=/){split(\$i,b,\"=\"); k=int(b[2])}}; print c, k}" $E/ptcl2D.txt > $E/class_cluster.txt; echo "distinct class ids: $(cut -d\" \" -f1 $E/class_cluster.txt | sort -un | wc -l), distinct cluster ids: $(cut -d\" \" -f2 $E/class_cluster.txt | sort -un | wc -l)"; cut -d" " -f1 $E/class_cluster.txt > $E/../subclass_labels.txt; done; echo DUMP DONE

#!/bin/bash
#SBATCH -J projcheck
#SBATCH -p csbdevel
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 00:20:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# integrity of a shared half project: jobproc entries, ptcl3D rows vs the pose table (args: PROJ POSETABLE)
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH; P=$1; T=$2; cd $(mktemp -d)
echo "== jobproc of $P"; simple_exec prg=print_project_field projfile=$P oritype=jobproc mkdir=no 2>&1 | grep -o "prg=[a-z0-9_]*" | sort | uniq -c
simple_exec prg=print_project_field projfile=$P oritype=ptcl3D mkdir=no 2>/dev/null > p3d.txt; echo "ptcl3D rows $(grep -c . p3d.txt), table rows $(grep -c . $T)"
for k in e1 e2 e3 x y dfx; do a=$(grep -o " $k=[^ ]*" p3d.txt | md5sum | cut -c1-12); b=$(grep -o " $k=[^ ]*" $T | md5sum | cut -c1-12); echo "$k project $a table $b $( [ $a = $b ] && echo same || echo DIFFERENT)"; done

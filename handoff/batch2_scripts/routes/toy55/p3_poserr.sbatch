#!/bin/bash
#SBATCH -J klem_p3poserr
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# pose errors of the rflex_s30a flex states vs truth (problem 3, p3_poserr.py)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K_/p3poserr; mkdir -p $O; S=$T/simswing30
for nd in r22 r2 r; do N=$K_/rflex_s30a/node_$nd; [ -s $O/dump_$nd.txt ] || simple_exec prg=print_project_field projfile=$N/n.simple oritype=ptcl3D mkdir=no 2>/dev/null > $O/dump_$nd.txt
  $PY $T/p3_poserr.py $O/dump_$nd.txt $N/flex_states.txt $S/all.txt --gt $S/gt_particles.txt --name s30a_$nd 2>&1 | grep -v Warning; done
echo DONE

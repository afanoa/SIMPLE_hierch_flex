#!/bin/bash
#SBATCH -J klem_landscape_cmp2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# port3's stage-B leaf labels (print_project_field), then landscape_cmp.sbatch's comparison with all four r1 runs + the r2 control
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH; K=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K/repro_r1/dscore_port3; mkdir -p $O
[ -s $O/labels.txt ] || simple_exec prg=print_project_field projfile=$K/stree_port/tree_f692r1_rep3/1_hetero_tree/node_r/stageB/b.simple oritype=ptcl3D mkdir=no 2>/dev/null \
  | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/labels.txt
echo "port3 labels: $(awk '$1 > 0' $O/labels.txt | wc -l) rows in $(awk '$1 > 0' $O/labels.txt | sort -u | wc -l) states"
bash /mnt/beegfs/elmlund/afan/routes/toy55/landscape_cmp.sbatch

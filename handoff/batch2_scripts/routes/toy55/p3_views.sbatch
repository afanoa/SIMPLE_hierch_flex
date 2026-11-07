#!/bin/bash
#SBATCH -J klem_p3views
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# view / defocus distributions of node states (problem 3, p3_views.py): real (rflex_f692rf flex states r, r2; f692r fam1 units) vs synthetic
# (rflex_s30a flex states r, r2; panel swing30hs refined states). Reads other runs' projects, writes only rnasep_klem/p3views.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K_/p3views; mkdir -p $O && cd $O
one(){ # name project labels
  [ -s $O/dump_$1.txt ] || simple_exec prg=print_project_field projfile=$2 oritype=ptcl3D mkdir=no 2>/dev/null > $O/dump_$1.txt
  $PY $T/p3_views.py $O/dump_$1.txt $3 --name $1 2>&1 | grep -v Warning; }
for nd in r r2; do one real_f692rf_$nd $K_/rflex_f692rf/node_$nd/n.simple $K_/rflex_f692rf/node_$nd/flex_states.txt; one sim_s30a_$nd $K_/rflex_s30a/node_$nd/n.simple $K_/rflex_s30a/node_$nd/flex_states.txt; done
one real_f692r_fam1 $K_/f692r/coordtree_fam1/rec/1_reconstruct3D/r.simple $K_/f692r/coordtree_fam1/units.txt
one sim_panel_swing30hs $K_/panel_swing30hs/1_refine3D_states/p.simple $K_/panel_swing30hs/refined_states.txt
echo DONE

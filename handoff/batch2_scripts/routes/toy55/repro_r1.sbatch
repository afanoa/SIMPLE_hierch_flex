#!/bin/bash
#SBATCH -J klem_repro_r1
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# run-to-run reproducibility of node partitions from existing runs (repro_r1.py): label files per run and step into rnasep_klem/repro_r1/labels
# (text label files copied, refine projects dumped with print_project_field), then the report rnasep_klem/repro_r1/report.txt
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; T=/mnt/beegfs/elmlund/afan/routes/toy55
O=$K_/repro_r1; L=$O/labels; mkdir -p $L; cd $O
S=$K_/rflex_f692st/node_r1; P=$K_/stree_port/tree_f692r1/1_hetero_tree/node_r; U=$K_/rflex_f692u7/node_r; N=$K_/rflex_f692u7n13f/node_r; BB=$K_/rflex_f692u7b13f/node_r
st(){ [ -s $L/$2.txt ] || simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $L/$2.txt; }
cp -f $S/flex_states.txt $L/S_flex.txt; cp -f $S/fflex_states.txt $L/S_fflex.txt; cp -f $S/refined_f.txt $L/S_B.txt; cp -f $S/old_cheap20/fflex_states.txt $L/O_fflex.txt
cp -f $U/flex_states.txt $L/U_flex.txt; cp -f $U/fflex_states.txt $L/U_fflex.txt; cp -f $U/refined_f.txt $L/U_B.txt
cp -f $N/flex_states.txt $L/N_flex.txt; cp -f $N/fflex_states.txt $L/N_fflex.txt; cp -f $N/refined_f.txt $L/N_B.txt; cp -f $BB/refined_f.txt $L/Bb_B.txt
st $S/fstageA/1_refine3D/r.simple S_A & st $S/fstageB/2_refine3D/b.simple S_Bdump & st $S/old_cheap20/fstageA/1_refine3D/r.simple O_A &
st $P/flex/f.simple P_flex & st $P/fflex/f.simple P_fflex & st $P/stageA/r.simple P_A & st $P/stageB/b.simple P_B & wait
st $U/fstageA/1_refine3D/r.simple U_A & st $N/fstageA/1_refine3D/r.simple N_A & st $BB/fstageA/1_refine3D/r.simple Bb_A & wait
for f in $L/*.txt; do printf "%-14s %8d lines, %8d labelled\n" $(basename $f .txt) $(wc -l < $f) $(awk '$1 > 0' $f | wc -l); done
cmp -s <(awk '{print ($1 > 0)}' $L/S_B.txt) <(awk '{print ($1 > 0)}' $L/S_Bdump.txt) && echo "refined_f.txt and the stage-B project agree on the node rows" || echo "WARNING: refined_f.txt and the stage-B project differ"
$PY $T/repro_r1.py $L > $O/report.txt 2>&1; cat $O/report.txt

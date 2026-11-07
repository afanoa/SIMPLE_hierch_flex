#!/bin/bash
#SBATCH -J klem_p3blur
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# local blur / scale regression before the focus mask (problem 3, p3_blur.py): realistic-SNR synthetic flex-state nodes with truth (rflex_s30a), the
# high-SNR panel root and two inner nodes' worth of leaves (all leaves = one node), against the real flex-state nodes and the f692r fam1 root (read only).
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; O=$K_/p3blur; mkdir -p $O && cd $O
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2; ARM=$T/maps160/toy_consensus.mrc:$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
mklist(){ # rdir labels nstates -> stdout
  awk -v d=$1 -v k=$3 'NR == FNR {c[$1]++; next} END {for (s = 1; s <= k; s++) printf "%d %s/recvol_state%02d_even_unfil.mrc %s/recvol_state%02d_odd_unfil.mrc %d\n", s, d, s, d, s, c[s]}' $2 /dev/null; }
P=$K_/panel_swing30hs; mklist $P/1_refine3D_states $P/refined_states.txt 20 > $O/list_panel_swing30hs.txt
C=$K_/f692r/coordtree_fam1; mklist $(ls -td $C/rec/*_reconstruct3D | head -1) $C/units.txt 30 > $O/list_f692r_fam1.txt
run(){ $PY $T/p3_blur.py $O "$@" > $O/log_$2.txt 2>&1 & }
for nd in r r2 r22 r12; do run --name s30a_$nd --list $K_/rflex_s30a/node_$nd/members.txt --labels $K_/rflex_s30a/node_$nd/flex_states.txt --gt $T/simswing30/gt_particles.txt --truth $T/simswing30/maps --mskdiam 220; done
run --name panel_swing30hs --list $O/list_panel_swing30hs.txt --labels $P/refined_states.txt --gt $T/simswing30hs/gt_particles.txt --truth $T/simswing30hs/maps --mskdiam 220
for nd in r r2; do run --name f692rf_$nd --list $K_/rflex_f692rf/node_$nd/members.txt --arm $ARM --mskdiam 232; done
run --name f692r_fam1 --list $O/list_f692r_fam1.txt --arm $ARM --mskdiam 232
wait; for f in $O/*_blur.txt; do cat $f; echo; done; grep -l "Traceback\|Error" $O/log_*.txt && echo "FAILURES above"; echo DONE

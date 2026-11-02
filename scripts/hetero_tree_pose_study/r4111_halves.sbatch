#!/bin/bash
#SBATCH -J r4111_halves
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/pose_study/logs/%x_%j.out
# leaf r4111_s1 of the 692k tree reconstructed from its top-corr half and its bottom-corr half at its stage-B poses (taken from
# the merge-up's master project, which carries every leaf's stage-B poses), both read with the saved readouts' chain
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | sed -n 's/^\(SIMPLE_[A-Z0-9_]*\)=.*/\1/p'); do unset $v; done
B=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/dist_sdisc/build12/install
export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH OMP_NUM_THREADS=16 PYTHONDONTWRITEBYTECODE=1
PY=/mnt/nasapps/production/python/3.12/bin/python3
H=/mnt/beegfs/elmlund/afan/hetero_tree_20261028; S=$H/fixtures/scripts; R160=/mnt/beegfs/elmlund/afan/routes/toy55/maps160
T=$H/integration/f692_r17d_s1/1_hetero_tree/node_r4111; M=$H/integration/f692_r17d_s1/1_hetero_tree/mergeup/final.simple
W=$H/pose_study/r4111_halves2
test -e $M || { echo "no $M"; exit 2; }
mkdir -p $W && cd $W
simple_exec prg=print_project_field projfile=$T/stageB/b.simple oritype=ptcl3D mkdir=no 2>/dev/null | grep 'e1=' > b_dump.txt
simple_exec prg=print_project_field projfile=$M oritype=ptcl3D mkdir=no 2>/dev/null | grep -c 'e1=' > nr.txt; echo "master rows $(cat nr.txt)"
simple_exec prg=print_project_field projfile=$M oritype=stk mkdir=no 2>/dev/null | head -1 | tr ' ' '\n' | grep -E '^(box|smpd)=' > stk.txt; cat stk.txt
$PY $W/../r4111_halves_flags.py $T $W
box=$(grep '^box=' stk.txt | cut -d= -f2 | cut -d. -f1); smpd=$(grep '^smpd=' stk.txt | cut -d= -f2)
smpdc=$(awk -v b=$box -v s=$smpd 'BEGIN{printf "%.5f", b*s/88.}')
echo "box $box smpd $smpd -> box_crop 88 smpd_crop $smpdc"
for h in A B; do
  mkdir -p $W/$h && cd $W/$h && cp $M l.simple
  simple_exec prg=selection projfile=l.simple oritype=ptcl3D infile=$W/sel$h.txt mkdir=no > selection.log 2>&1
  echo "-- $h: selection errors $(grep -c ERROR selection.log); state-1 rows after selection: $(simple_exec prg=print_project_field projfile=l.simple oritype=ptcl3D mkdir=no 2>/dev/null | grep -c ' state=1')"
  simple_exec prg=reconstruct3D projfile=l.simple pgrp=c1 mskdiam=232 nthr=16 box_crop=88 smpd_crop=$smpdc ml_reg=no postprocess=no mkdir=no > reconstruct3D.log 2>&1
  echo "-- $h: $(grep -c 'NORMAL STOP' reconstruct3D.log) normal stop; skipped: $(grep -h 'SKIPPING' reconstruct3D.log | tail -1); maps: $(find $W/$h -name 'recvol_state01*_unfil.mrc' | wc -l)"
done
cd $W
E_A=$(find $W/A -name 'recvol_state01_even_unfil.mrc' | head -1); O_A=$(find $W/A -name 'recvol_state01_odd_unfil.mrc' | head -1)
E_B=$(find $W/B -name 'recvol_state01_even_unfil.mrc' | head -1); O_B=$(find $W/B -name 'recvol_state01_odd_unfil.mrc' | head -1)
printf '%s\n' "r4111_s1_top $E_A $O_A" "r4111_s1_bottom $E_B $O_B" "r4111_s1_full $T/stageB/recvol_state01_even_unfil.mrc $T/stageB/recvol_state01_odd_unfil.mrc" > spec.txt
$PY $S/make_view.py $W/views $W/spec.txt --box 160 2>&1 | grep -v Deprecation | tail -2
CL=$(ls $R160/docked/conf18_*_m160.mrc | head -1)
$PY $S/arm_readout.py $R160/toy_consensus.mrc $CL $W/arm.txt $W/views/r4111_s1_top.mrc $W/views/r4111_s1_bottom.mrc $W/views/r4111_s1_full.mrc 2>&1 | grep -v Deprecation | tail -2
echo "== arm.txt"; cat $W/arm.txt
md5sum $E_A $E_B

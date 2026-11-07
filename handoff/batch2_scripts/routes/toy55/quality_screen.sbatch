#!/bin/bash
#SBATCH -J quality_screen
#SBATCH -p norm
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/scorecard/%x.%j.out
# Validation of the knob-free quality screen (quality_screen.py) on cases with known answers, existing maps only:
# - quality: the m33 root's site 6 (gate/m33_root_m79, bins 1 vs 5); the 22 Oct tree's r3 and r4 against r1 (node means);
# - real: root sites 1 and 2; the 22 Oct r2 vs r1 and its largest closed-arm leaf vs its largest open r2 leaf; cdisc closed and A,
#   IgM a4_th and a2_th, CAHRA open (phase E bins 1 vs 5); unlabelled: root sites 3, 5 and 7.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_stree/snap/m8c52t; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local SIMPLE_EMAIL=nil OMP_NUM_THREADS=1
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port; E=$K/factorial/phaseE
G=$S/gate/m33_root_m79/1_hetero_tree/node_r; H22=$S/full692/1_hetero_tree; F10=$K/full692_analysis/final10
CD=$S/fulltoy_m52M_cd1/1_hetero_tree; IG=$S/fulltoy_m52M_igm/1_hetero_tree; CA=$S/arm_M_c11_s1/1_hetero_tree
OUT=$K/scorecard/quality_screen; mkdir -p $OUT/g; cd $OUT
dump(){ [ -s $2 ] || simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D 2>/dev/null | $PY $T/dump_ptcl.py $2 > $2.log 2>&1; }
dump $G/stageB/b.simple d_m79.npz & dump $H22/node_r/stageB/b.simple d_f692.npz & dump $CD/node_r/stageB/b.simple d_cd.npz &
dump $IG/node_r/stageB/b.simple d_igm.npz & dump $CA/node_r/stageB/b.simple d_cahra.npz & wait
nr(){ $PY -c "import numpy as np, sys; print(len(np.load(sys.argv[1])['corr']))" $1; }
mg(){ $PY $T/make_groups.py "$@" >> groups.log 2>&1; }
C=cases.txt; : > $C; N79=$(nr d_m79.npz); N22=$(nr d_f692.npz)
for c in c1 c2 c3 c5 c6 c7; do M=$G/motion/$c; case $c in c6) tr=quality;; c1|c2) tr=real;; *) tr=unknown;; esac
  mg $N79 g/m79_$c.txt bins:$M/bins.txt:1 bins:$M/bins.txt:5
  echo "m33root_site${c#c} $tr $M/recvol_state01.mrc $M/recvol_state01_even.mrc $M/recvol_state01_odd.mrc $M/recvol_state05.mrc $M/recvol_state05_even.mrc $M/recvol_state05_odd.mrc $G/mask.mrc d_m79.npz g/m79_$c.txt" >> $C; done
nm(){ echo "$H22/node_$1/mean/recvol_state01.mrc $H22/node_$1/mean/recvol_state01_even_unfil.mrc $H22/node_$1/mean/recvol_state01_odd_unfil.mrc"; }
for x in r3 r4 r2; do tr=quality; [ $x = r2 ] && tr=real; mg $N22 g/f692_${x}_r1.txt rows:$H22/node_r1/rows.txt rows:$H22/node_$x/rows.txt
  echo "oct22_${x}_vs_r1 $tr $(nm r1) $(nm $x) $H22/node_r/mask.mrc d_f692.npz g/f692_${x}_r1.txt" >> $C; done
lc=$(awk '!/^#/ && $2 == "leaf" && $5 >= 0.6 {print $4, $1}' $F10/maps_table.txt | sort -n | tail -1 | cut -d' ' -f2)
lo=$(awk '!/^#/ && $2 == "leaf" && $3 == "r2" && $5 <= 0.3 {print $4, $1}' $F10/maps_table.txt | sort -n | tail -1 | cut -d' ' -f2)
kc=$(awk -v n=$lc '$1 == n {print NR}' $H22/hetero_tree_leaves.txt); ko=$(awk -v n=$lo '$1 == n {print NR}' $H22/hetero_tree_leaves.txt)
mg $N22 g/f692_closed_open.txt eq:$H22/hetero_tree_labels.txt:$ko eq:$H22/hetero_tree_labels.txt:$kc
echo "oct22_closed_${lc}_vs_open_${lo} real $(awk -v n=$lo '$1 == n {print $5, $6, $7}' $F10/maplist.txt) $(awk -v n=$lc '$1 == n {print $5, $6, $7}' $F10/maplist.txt) $H22/node_r/mask.mrc d_f692.npz g/f692_closed_open.txt" >> $C
ser(){ local nm=$1 d=$2 mk=$3 dp=$4; mg $(nr $dp) g/$nm.txt bins:$d/work/states.txt:1 bins:$d/work/states.txt:5
  echo "$nm real $d/bin1.mrc $d/halves/bin1_even.mrc $d/halves/bin1_odd.mrc $d/bin5.mrc $d/halves/bin5_even.mrc $d/halves/bin5_odd.mrc $mk $dp g/$nm.txt" >> $C; }
ser cdisc_closed $E/switch2/cdisc_closed_truthjunk $CD/node_r/mask.mrc d_cd.npz; ser cdisc_A $E/switch/cdisc_A_truthjunk $CD/node_r/mask.mrc d_cd.npz
ser igm_a4th $E/igm/m52M_r_s1_p1_a4th $IG/node_r/mask.mrc d_igm.npz; ser igm_a2th $E/igm/m52M_r2_s4_p1_a2th $IG/node_r/mask.mrc d_igm.npz
ser cahra_open $E/switch/cahra_open_bc $CA/node_r/mask.mrc d_cahra.npz
echo "[$(date "+%F %T") quality_screen] $(wc -l < $C) cases"
split -n l/8 -d $C part_
for f in part_0*; do [ -s $f ] && ( $PY $T/quality_screen.py $f > $f.out 2>&1 ) & done; wait
cat part_0*.out | grep "^QSCREEN" > qscreen.txt; grep -h -i "error\|traceback" part_0*.out | head -5
echo "[$(date "+%F %T") quality_screen] DONE: $(wc -l < qscreen.txt) of $(wc -l < $C) cases"

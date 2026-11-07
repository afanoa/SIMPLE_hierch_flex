#!/bin/bash
#SBATCH -J klem_calF_rel
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 00:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# phaseF_rel.py (how big and how local each series' end-to-end change is) on every phase-E series (calibration, nulls, RNase P m33),
# the phase-B / C series they rescore and OLE's three, 8 lists in parallel; RELF lines in factorial/phaseE/relF.txt; phaseE on OLE too
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; T=/mnt/beegfs/elmlund/afan/routes/toy55; K=/mnt/beegfs/elmlund/afan/rnasep_klem; S=$K/stree_port; E=$K/factorial/phaseE
IG=$S/fulltoy_m52M_igm/1_hetero_tree; MB=$S/fulltoy_igm/1_hetero_tree; CA=$S/arm_M_c11_s1/1_hetero_tree; CD=$S/fulltoy_m52M_cd1/1_hetero_tree; M33=$S/full692_m33_s1/1_hetero_tree
I=$K/factorial/motion_series_igm; PB=$K/factorial/phaseB; O=/mnt/beegfs/elmlund/afan/ole_12707; L=()
ok(){ [ -s $1/bin5.mrc ] && [ -s $1/halves/bin5_odd.mrc ]; }
for d in $E/igm/m52M_*; do ok $d && L+=("$(basename $d):$d:$IG/node_r/mask.mrc"); done
for nm in m52M_r_s1_a4th_q5 m52M_r2_s4_a2th_q5; do L+=("ex_$nm:$I/$nm:$IG/node_r/mask.mrc"); done; L+=("ex_m8c17_r2_s1_a2th_q5:$I/m8c17_r2_s1_a2th_q5:$MB/node_r/mask.mrc")
while read _ ang site j rho f lev nr x; do [ "$x" = same ] && continue; nm=m8c17_${site}_p${j}_${ang/_/}; ok $E/igm_m8c17/$nm && L+=("$nm:$E/igm_m8c17/$nm:$MB/node_$lev/mask.mrc"); done < <(grep "^PICK" $E/igm_m8c17/picks2.txt)
ok $E/switch/cahra_open_bc && L+=("cahra_open_bc:$E/switch/cahra_open_bc:$CA/node_r/mask.mrc"); L+=("ex_cahra_open:$PB/cahra_open:$CA/node_r/mask.mrc")
for d in $E/switch/cdisc_* $E/switch2/cdisc_*; do ok $d && L+=("$(basename $d):$d:$CD/node_r/mask.mrc"); done
for nm in cdisc_A cdisc_closed; do L+=("ex_$nm:$PB/$nm:$CD/node_r/mask.mrc"); done
for d in $E/null/null_* $E/null2/null_* $E/null3/null_*; do ok $d || continue; b=$(basename $d); nd=${b##*_}; case $b in null_igm_*) L+=("$b:$d:$IG/node_$nd/mask.mrc");; null_cdisc_*) L+=("$b:$d:$CD/node_$nd/mask.mrc");; esac; done
for d in $E/rnasep_m33/node_*/m33_*; do ok $d || continue; nd=$(basename $(dirname $d)); L+=("$(basename $d):$d:$M33/$nd/mask.mrc"); done
for s in r1_s4 r2_s4 r_s2; do L+=("ole_$s:$O/motion_series/$s:$O/6_hetero_tree/node_r/mask.mrc"); done
echo "[$(date "+%F %T") calF_rel] ${#L[@]} series"; mkdir -p $E/relF_parts
( OMP_NUM_THREADS=1 $PY $T/phaseE_maptests.py $(for s in r1_s4 r2_s4 r_s2; do printf "ole_%s:%s/motion_series/%s:%s/6_hetero_tree/node_r/mask.mrc:240 " $s $O $s $O; done) > $E/ole_maptests.txt 2>&1 ) &
for p in 0 1 2 3 4 5 6 7; do A=(); for ((i = p; i < ${#L[@]}; i += 8)); do A+=("${L[$i]}"); done
  [ ${#A[@]} -gt 0 ] && ( OMP_NUM_THREADS=1 $PY $T/phaseF_rel.py "${A[@]}" > $E/relF_parts/part$p.txt 2>&1 ) &
done; wait
cat $E/relF_parts/part*.txt | grep -v Warning > $E/relF.txt; echo "[$(date "+%F %T") calF_rel] DONE: $(grep -c ^RELF $E/relF.txt) RELF lines, $(grep -c -i -E "error|traceback" $E/relF.txt) error lines"

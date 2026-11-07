#!/bin/bash
#SBATCH -J klem_p1hold
#SBATCH -p csbdevel
#SBATCH -c 48
#SBATCH --mem=200G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# held-out comparison of a coordinate node's leaf models (problem 1): the node's hard states, NB equal-width bins of its motion
# coordinate and a random split, each reconstructed at the node's poses; every particle scored against the OTHER half's FSC-weighted
# maps of all 3*NB partitions (refine3D maxits=1 with lp set matches both halves against volN, so two passes over all particles:
# pass E keeps the even rows scored against odd-half maps, pass O the odd rows); SELF=1 adds pass S = the hard maps as averages
# (self-inclusive, as the coordinate stage scored). Analysis p1_holdout.py. env ND (node), NB (8), RT (s30a), TAG (swing30), LP (8), SELF (0)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=12 OPENBLAS_NUM_THREADS=12
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
ND=${ND:?}; NB=${NB:-8}; RT=${RT:-s30a}; TAG=${TAG:-swing30}; GT=$T/sim$TAG/gt_particles.txt; R=$K_/rflex_$RT; C=$R/coord_$ND
O=$K_/p1hold/$ND; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') p1hold $ND] $*"; }
[ -s base_dump.txt ] || { cp -f $C/base.simple base.simple && simple_exec prg=print_project_field projfile=base.simple oritype=ptcl3D mkdir=no 2>/dev/null > base_dump.txt; }
$PY $T/p1_parts.py $O --hard $R/node_$ND/refined_states.txt --coord $C/motion_coord.txt --dump base_dump.txt --nb $NB || { log "partitions FAILED"; exit 1; }
rec(){ # $1 partition -> $O/rec_$1/1_reconstruct3D
  mkdir -p $O/rec_$1 && cd $O/rec_$1; cp -f $O/base.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/$1.txt mkdir=no > rsel.log 2>&1 || { log "selection $1 FAILED"; exit 1; }
  simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=12 mkdir=no > rup.log 2>&1
  simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=220 nstates=$NB nparts=4 nthr=12 > rec.log 2>&1 || { log "reconstruct3D $1 FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; exit 1; }
  RR=$(ls -td $O/rec_$1/*_reconstruct3D | head -1); rm -f $RR/*_part*_even.mrc $RR/*_part*_odd.mrc $RR/rho_*; cd $O; }
for m in hard bins rand; do [ -n "$(ls $O/rec_$m/*_reconstruct3D/recvol_state$(printf %02d $NB)_odd_unfil.mrc 2>/dev/null)" ] || { log "reconstructing $m"; rec $m; }; done
$PY $T/p1_refs.py $O $(for m in hard bins rand; do echo --rec $(ls -td $O/rec_$m/*_reconstruct3D | head -1); done) --nb $NB --self $([ "${SELF:-0}" = 1 ] && echo $NB || echo 0) || { log "references FAILED"; exit 1; }
score(){ # $1 pass, $2 state file, $3 nstates, $4 refs dir
  W=$O/pass_$1; [ -n "$(ls $W/*_refine3D/dist_part_neigh_4.dat 2>/dev/null)" ] && return 0; rm -rf $W; mkdir -p $W && cd $W; cp -f $O/base.simple s.simple
  simple_exec prg=selection projfile=s.simple oritype=ptcl3D infile=$2 mkdir=no > sel.log 2>&1 || { log "pass $1 selection FAILED"; exit 1; }
  simple_exec prg=update_project projfile=s.simple qsys_name=local nparts=4 nthr=12 mkdir=no > up.log 2>&1
  V=""; for j in $(seq 1 $3); do V="$V vol$j=$4/vol_$(printf %02d $j).mrc"; done
  log "pass $1: $3 maps from $(basename $4), lp ${LP:-8}"
  simple_exec prg=refine3D projfile=s.simple $V nstates=$3 refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=220 lp=${LP:-8} maxits=1 center=no nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=12 > score.log 2>&1 || { log "pass $1 refine3D FAILED: $(grep -a -m2 ERROR score.log | cut -c1-150)"; exit 1; }
  SD=$(ls -td $W/*_refine3D | head -1); rm -f $SD/*_part*_even.mrc $SD/*_part*_odd.mrc $SD/rho_* $SD/*.mrc; cd $O; }
score E $O/score_states.txt $((3 * NB)) $O/refsE; score O $O/score_states.txt $((3 * NB)) $O/refsO
[ "${SELF:-0}" = 1 ] && score S $O/score_states_hard.txt $NB $O/refsS
$PY $T/p1_holdout.py $O --passE $O/pass_E --passO $O/pass_O $([ "${SELF:-0}" = 1 ] && echo --passS $O/pass_S) --nb $NB $([ -n "$TAG" ] && echo --gt $GT) --it0 $C/it0/pc 2>&1 | grep -v Warning
log DONE

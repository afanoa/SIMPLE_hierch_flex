#!/bin/bash
#SBATCH -J klem_focused_arm
#SBATCH -p norm
#SBATCH -c 64
#SBATCH --mem=240G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# arm-focused classification of one real rflex node: (1) real-space variation mask (arm_mask.py) from the node's view-balanced children along its
# arm axis (cheapkids rec_1: 10 quantile bins of PC1); (2) refine3D of the node's particles at the parent's poses, seeded by the EQUAL-WIDTH
# bins along the same axis (rec_w1: 12 bins, small end bins), pcg_focus with that mask (node mean fixed outside), regular settings;
# (3) closed-arm readout per state (view maps lp 10 in view_rflex/focused_arm_<RT>_<NODE>/) and where an earlier run's closed states land.
# env RT (f692rf), NODE (r22), ITS (6), LP (10), MSKDIAM (232), FRAC (0.2)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:-f692rf}; NODE=${NODE:-r22}; MD=${MSKDIAM:-232}
R=$K_/rflex_$RT; CK=$K_/cheapkids_${RT}_$NODE; O=$K_/focused_arm_${RT}_$NODE; mkdir -p $O && cd $O; log(){ echo "[$(date '+%F %T') focused-arm $RT $NODE] $*"; }
CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); set -- $(awk -v u=$NODE '$1 == u' $R/queue.txt); PP=${3:-}
R1=$(ls -td $CK/rec_1/*_reconstruct3D | head -1); RW=$(ls -td $CK/rec_w1/*_reconstruct3D | head -1); NW=$(awk '$1 > 0' $CK/children_w1.txt | sort -n | uniq | wc -l)
[ -s "$PP" ] && [ -d "$R1" ] && [ -d "$RW" ] || { log "inputs missing"; exit 1; }
: > $O/mask_list.txt; for s in $(seq 1 10); do b=$R1/recvol_state$(printf %02d $s); echo "c$s ${b}_even_unfil.mrc ${b}_odd_unfil.mrc $(awk -v s=$s '$1 == s' $CK/children_1.txt | wc -l)" >> $O/mask_list.txt; done
[ -s $O/arm_mask.mrc ] || $PY $T/arm_mask.py $O/arm_mask.mrc --list $O/mask_list.txt --mskdiam $MD --frac ${FRAC:-0.2} --cons $T/maps160/toy_consensus.mrc --closed $CL 2>&1 | grep -v Warning
VOLS=""; for s in $(seq 1 $NW); do VOLS="$VOLS vol$s=$RW/recvol_state$(printf %02d $s).mrc"; done
if [ -z "$(ls $O/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cp -f $PP r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$CK/children_w1.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=local nparts=4 nthr=16 mkdir=no > up.log 2>&1 || { log "project FAILED"; exit 1; }
  log "arm-focused refine3D: $NW states seeded by the equal-width arm-axis bins ($(awk '$1 > 0' $CK/children_w1.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')), lp ${LP:-10}, maxits ${ITS:-6}"
  simple_exec prg=refine3D projfile=r.simple $VOLS nstates=$NW refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=${LP:-10} maxits=${ITS:-6} center=no rec_backend=pcg pcg_mskfile=$O/arm_mask.mrc pcg_focus=yes nspace=5000 box_crop=88 smpd_crop=2.66181803 nparts=4 nthr=16 > ref.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }; fi
RD=$(ls -td $O/*_refine3D | head -1); [ -s $O/states.txt ] || simple_exec prg=print_project_field projfile=$RD/r.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}' > $O/states.txt
log "state populations: $(awk '$1 > 0' $O/states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}'); overlap per iteration: $(grep -a 'STATE OVERLAP' $RD/../ref.log 2>/dev/null | awk '{printf "%s ", $NF}')"
V=$K_/view_rflex/focused_arm_${RT}_$NODE; : > $O/spec.txt; for s in $(seq 1 $NW); do b=$RD/recvol_state$(printf %02d $s); echo "state_$(printf %02d $s) ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec.txt; done
$PY $T/make_view.py $V $O/spec.txt --box 160 > /dev/null 2>&1
$PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/arm.txt $(ls $V/*.mrc) 2>&1 | grep -v Warning | tail -1
$PY - <<PYEOF
import numpy as np, re
st = np.loadtxt('$O/states.txt', dtype=int); ref = np.loadtxt('$K_/f692r/coordtree_fam1/units.txt', dtype=int)
arm = {int(re.findall(r'state_(\d+)', l.split()[0])[0]): (float(l.split()[1]), float(l.split()[2])) for l in open('$O/arm.txt') if l.strip() and not l.startswith('#')}
print('state: particles, closed-arm, open-arm | share of the earlier closed states 16 / 15 (closed-arm 0.87 / 0.83) in this state, and the state\'s purity in them')
for s in sorted(arm, key=lambda s: -arm[s][0]):
    m = st == s; c16 = (ref[m] == 16).sum(); c15 = (ref[m] == 15).sum()
    print('   %2d: %6d  %.2f %.2f | %4.0f %% / %4.0f %% of states 16 / 15; %.1f %% of the state' % (s, m.sum(), arm[s][0], arm[s][1], 100 * c16 / max((ref[st > 0] == 16).sum(), 1), 100 * c15 / max((ref[st > 0] == 15).sum(), 1), 100 * (c16 + c15) / max(m.sum(), 1)))
PYEOF
log DONE

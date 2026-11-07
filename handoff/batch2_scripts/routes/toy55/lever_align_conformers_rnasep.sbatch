#!/bin/bash
#SBATCH -J rnasep_lever_conformers
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The settled node policy applied to every identified conformer of family 1: for each of the ten matched states of the
# feedback board (matched_rec/states.txt, state s = conformer ids[s]), a one-cell klem board of its particles with the
# whole-reference local pose step (rotation-only, 3x3x3 grid x 2 refinements, Gaussian prior), the moved poses of all
# ten sets imported into one fresh project (stack order, ctf=flip, explicit selection with the state labels) and
# reconstructed as ten states; identification against the deposited conformers before/after and the FSC per map.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; P=$E/fam1_pipeline; KD=$P/klem_fb; META=$P/meta_fb.npz
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
STATES=$KD/matched_rec/states.txt; D=$KD/lever_conformers; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') lc] $*"; }; IDS="1 10 18 27 28 37 46 48 49 51"
NS=$(awk '$1>0{print $1}' $STATES | sort -n | uniq | wc -l); log "$NS matched states (conformers $IDS)"
export OMP_NUM_THREADS=1
for s in $(seq 1 $NS); do
  mkdir -p $D/s$s; awk -v s=$s '$1==s{print NR}' $STATES > $D/s$s/subset.txt; n=$(grep -c . $D/s$s/subset.txt); [ $n -ge 1000 ] || { log "state $s: $n particles, skipped"; continue; }
  [ -s $D/s$s/board/ptcl3D_pose_it6.txt ] && { log "state $s: poses exist"; continue; }
  log "state $s ($n particles): one-cell board with the whole-reference pose step"
  $PY $W/klem.py $STK none none $D/s$s/board --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 1 --n2 1 --nfree 0 --nproc 16 --init random --nsub 0 --subset $D/s$s/subset.txt --ctf auto \
      --lps "10 10 10 8 8 8" --sigmas "0.7 0.7 0.7 0.7 0.7 0.7" --mskdiam 232 --pose_every 2 --pose_step 3 --pose_prior 4 --pose_ref cell --tag "lc s$s" 2>&1 | grep "pose step (\|Error\|Traceback" | tail -2 | cut -c1-160
done
export OMP_NUM_THREADS=8; ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort > $D/stktab.txt
$PY - <<PYEOF
import numpy as np, glob
M = np.load('$META'); N = len(M['flags']); e1, e2, e3, x, y = (M[k].astype(float).copy() for k in ('e1', 'e2', 'e3', 'x', 'y')); dfx, dfy, ang = (M[k].astype(float) for k in ('dfx', 'dfy', 'angast')); stk = M['stkind'].astype(int); ind = M['indstk'].astype(int); eo = M['eo'].astype(int) if 'eo' in M else np.arange(N) % 2
state = np.zeros(N, int); moved = 0; d0 = []
for s in range(1, $NS + 1):
    fs = sorted(glob.glob('$D/s%d/board/ptcl3D_pose_it*.txt' % s), key=lambda f: int(f.split('_it')[-1].split('.')[0]))
    if not fs: continue
    rows = [dict(t.split('=', 1) for t in l.split()) for l in open(fs[-1]) if l.strip()]
    for r in rows:
        p = int(r['pind']) - 1; d0.append(abs(float(r['e2']) - e2[p])); e1[p], e2[p], e3[p] = float(r['e1']), float(r['e2']), float(r['e3']); state[p] = s; moved += 1
order = np.lexsort((ind, stk)); n_stk = [int(np.fromfile(f.strip(), np.int32, 256)[2]) for f in open('$D/stktab.txt')]
print('%d particles moved over %d states (|d e2| median %.2f deg); labels %s; stacks %s' % (moved, $NS, np.median(d0) if d0 else 0., np.bincount(state)[1:], n_stk))
with open('$D/oritab.txt', 'w') as f:
    for i in order: f.write('e1=%.4f e2=%.4f e3=%.4f x=%.3f y=%.3f dfx=%.5f dfy=%.5f angast=%.3f eo=%d state=%d\n' % (e1[i], e2[i], e3[i], x[i], y[i], dfx[i], dfy[i], ang[i], eo[i], state[i]))
np.savetxt('$D/states_stackorder.txt', state[order], fmt='%d')
PYEOF
cd $D && rm -rf proj && simple_exec prg=new_project projname=proj > np.log 2>&1 && cd proj || { log "new_project FAILED"; exit 1; }
simple_exec prg=import_particles projfile=proj.simple stktab=$D/stktab.txt oritab=$D/oritab.txt ctf=flip smpd=1.464 kv=300 cs=2.7 fraca=0.1 mkdir=no > import.log 2>&1 || { log "import FAILED"; tail -3 import.log; exit 1; }
simple_exec prg=selection projfile=proj.simple oritype=ptcl3D infile=$D/states_stackorder.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
log "import + selection done: $(simple_exec prg=print_project_field projfile=proj.simple oritype=ptcl3D mkdir=no 2>/dev/null | grep -o ' state=[1-9][0-9]*\.' | sort | uniq -c | tr '\n' ' ')"
simple_exec prg=reconstruct3D projfile=proj.simple pgrp=c1 mskdiam=232 nstates=$NS nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-160)"; exit 1; }
RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $D/maps $D/states_dir; for s in $(seq -w 1 $NS); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $D/maps/ 2>/dev/null; done; cp -f $RR/recvol_state$s.mrc $D/states_dir/state_$s.mrc; done; cd $D
mkdir -p $D/states_dir_before; for s in $(seq -w 1 $NS); do cp -f $KD/matched_rec/maps/recvol_state$s.mrc $D/states_dir_before/state_$s.mrc 2>/dev/null; done
echo "== identification BEFORE (family poses, matched_rec maps)"; $PY $W/deposited_match.py $D/states_dir_before $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -13
echo "== identification AFTER (whole-reference pose step per state)"; $PY $W/deposited_match.py $D/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -13
VIS=$KD/closed_focal/lever/vis.py; args=(); for s in $(seq 1 $NS); do c=$(echo $IDS | cut -d" " -f$s); args+=("$KD/matched_rec/maps/recvol_state$(printf %02d $s).mrc" "conf $c before"); args+=("$D/maps/recvol_state$(printf %02d $s).mrc" "conf $c AFTER"); done; $PY $VIS "${args[@]}" 2>&1 | grep -v Warning
mkdir -p $D/for_download; for s in $(seq 1 $NS); do c=$(echo $IDS | cut -d" " -f$s); SRC=$D/maps/recvol_state$(printf %02d $s).mrc DST=$D/for_download/conf${c}_leveraligned_lp10_unsharpened.mrc BFAC=0 LPFIX=10 bash $W/sharpen_one_rnasep.sbatch > /dev/null 2>&1; done; ls $D/for_download | wc -l
log DONE

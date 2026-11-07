#!/bin/bash
#SBATCH -J rnasep_region_board
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Truth-free node recipe for ONE data-derived region: (1) focal klem board over the family with the references free
# only inside MASKFILE (region_k.mrc from derive_regions.py), at the family poses; (2) states = clusters of the board's
# reproducible cells (node_states.py), each state's particles ranked by their responsibility mass in the cluster, the
# top N reconstructed unmasked at the family poses; (3) level 2 for the largest states: a board of the state's own
# particles with klem's data-derived focus, cells ranked against the parent state's deviation (rank_by_target.py), the
# top N2 reconstructed. The deposited maps enter only in the identification readouts (deposited_match.py, vis.py).
# env: MASKFILE, REG (tag), SUBSET (family1_subset.txt), PIPE (fam1_pipeline), NMAX (30000), N2 (10000), NL2 (3), SIM (0.7), BAL (0.5)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; E=$A/rnasep_klem/expand1; P=${PIPE:-$E/fam1_pipeline}; KD=$P/klem_fb; META=$P/meta_fb.npz
CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt; STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
SUBSET=${SUBSET:-$E/family1_subset.txt}; D=$A/rnasep_klem/knobfree/region_${REG}; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') rb$REG] $*"; }; VIS=$E/fam1_pipeline/klem_fb/closed_focal/lever/vis.py
NROWS=$($PY -c "import numpy as np; print(len(np.load('$META')['flags']))"); PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple
board(){ out=$1; sub=$2; tag=$3; shift 3; export OMP_NUM_THREADS=1
  $PY $W/klem.py $STK none none $out --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 0 --nproc 16 --init npz --nsub 0 --subset $sub --ctf auto \
    --lps "12 12 10 10 8 8 8 8 8" --sigmas "1.2 1.2 1.0 1.0 0.8 0.8 0.7 0.7 0.7" --pose_marg 3 --devw 0.05 --pnoise --balance ${BAL:-0.5} --tag "$tag" "$@" 2>&1 | grep --line-buffered "subset\|states: cross\|focus mask\|Error\|Traceback" | cut -c1-200
  export OMP_NUM_THREADS=8; }
rec(){ dir=$1; sf=$2; ns=$3; mkdir -p $dir/rec && cd $dir/rec && cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$sf mkdir=no > sel.log 2>&1 || { log "selection FAILED ($dir)"; cd $D; return 1; }
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED ($dir): $(grep -a -m2 ERROR rec.log | cut -c1-200)"; cd $D; return 1; }
  RR=$(ls -td ./*_reconstruct3D | head -1); mkdir -p $dir/maps $dir/states_dir; for s in $(seq -f "%02g" 1 $ns); do for h in "" _even _odd; do cp -f $RR/recvol_state$s$h.mrc $dir/maps/ 2>/dev/null; done; cp -f $RR/recvol_state$s.mrc $dir/states_dir/state_$s.mrc; done; cd $D; }
log "(1) focal board on $(grep -c . $SUBSET) particles, mask $MASKFILE"; board $D/board $SUBSET "region $REG" --focus_mask $MASKFILE
log "(2) states from the board's reproducible cell clusters"; $PY $W/node_states.py $D/board $D/level1 --nrows $NROWS --sim ${SIM:-0.7} --nmax ${NMAX:-30000} 2>&1 | grep -v Warning
K=$(grep -vc "^#" $D/level1/clusters.txt); [ "$K" -ge 1 ] || { log "no states"; exit 1; }
rec $D/level1 $D/level1/states_clusters.txt $K || exit 1
echo "== identification of the level-1 states of region $REG (deposited maps used only here)"; $PY $W/deposited_match.py $D/level1/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc --save $D/level1/cell_matches.txt 2>&1 | grep -v Warning | tail -14
args=(); for s in $(seq 1 $K); do n=$(awk -v s=$s '$1==s' $D/level1/states_clusters.txt | wc -l); args+=("$D/level1/maps/recvol_state$(printf %02d $s).mrc" "region $REG state $s ($n ptcls)"); done; $PY $VIS "${args[@]}"
log "(3) level 2 for the largest ${NL2:-3} states"
for s in $(seq 1 ${NL2:-3}); do [ $s -le $K ] || break; L=$D/level2_s$s; mkdir -p $L; awk -v s=$s '$1==s{print NR}' $D/level1/states_clusters.txt > $L/subset.txt; n=$(grep -c . $L/subset.txt); [ $n -ge 4000 ] || { log "state $s: only $n particles, no level 2"; continue; }
  log "state $s: within-state board on $n particles, data-derived focus"; board $L/board $L/subset.txt "region $REG level2 s$s" --focus 1
  n2=$(( n * 4 / 10 )); [ $n2 -le ${N2:-10000} ] || n2=${N2:-10000}; $PY $W/rank_by_target.py $L/board $D/level1/cluster_targets.npz --target $s --n $n2 --nrows $NROWS --out $L 2>&1 | grep -v Warning
  rec $L $L/states_top.txt 1 || continue
  echo "== identification of region $REG state $s level 2 (top $n2)"; $PY $W/deposited_match.py $L/states_dir $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -12
  $PY $VIS $D/level1/maps/recvol_state$(printf %02d $s).mrc "region $REG state $s level 1 ($n)" $L/maps/recvol_state01.mrc "region $REG state $s level 2 (top $n2)"
done
log DONE

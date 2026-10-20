#!/bin/bash
#SBATCH -J toy10_family_feedback
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Lever 5: family-level pose feedback. For each family, the poses refined by the family's 10-state refinement go back
# into a child klem run on the family's particles (no free comps, best configuration); coverage per family and merged.
# env: NFAM (2), E (expansion dir with family subsets), RTAG of the refinement dirs
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; NFAM=${NFAM:-2}; E=${E:-$W/lat/expand2b}; CH=""
for j in $(seq 1 $NFAM); do
  RD=$(ls -td $W/lat/state_refine_k10_fam$j/*_refine3D_states | head -1); T=$W/lat/fb_fam$j; mkdir -p $T
  simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $T/refined_ptcl3D.txt
  awk '{ gsub(/state=[0-9.]+/, "state=1"); print }' $T/refined_ptcl3D.txt > $T/refined_state1.txt
  echo "[$(date '+%F %T')] family $j: poses from $RD ($(wc -l < $T/refined_state1.txt) rows); $($PY $W/pose_error.py $T/refined_ptcl3D.txt $W/sim/all.txt --tag poses)"
  env KTAG=_fb_fam$j INIT=pca NSUB=0 SUBSET=$E/family${j}_subset.txt NFREE=0 TABLE=$T/refined_state1.txt PMARG=3 DEVW=0.05 PNOISE=1 bash $W/run_klem.sbatch 2>&1 | grep "subset\|final:\|states: cross-half\|Error\|Traceback" | sed "s/^/[fam $j] /"
  CH="$CH $W/lat/klem_fb_fam$j"
done
$PY $W/expand_families.py merge $E $CH --gt $W/sim/gt_particles.txt --minrep 0.5
echo "[$(date '+%F %T')] feedback DONE"

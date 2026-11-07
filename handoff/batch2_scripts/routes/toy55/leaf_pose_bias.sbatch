#!/bin/bash
#SBATCH -J klem_posebias
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 00:40:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# per-state common pose error of a tree's merged states (env M = set tag e.g. null_50k_A, ST = states file (default
# consolidated_xtree/states.txt)), with the families' refine3D_states poses and with the initial consensus poses
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; M=${M:?}; R=$A/rnasep_klem/synth_$M; ST=${ST:-$R/consolidated_xtree/states.txt}
TABS=""; SUBS=""
for f in 1 2; do RD=$(ls -td $R/fam$f/refine/*_refine3D_states 2>/dev/null | head -1); [ -n "$RD" ] || continue
  [ -s $R/fam$f/refined_ptcl3D_print.txt ] || simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $R/fam$f/refined_ptcl3D_print.txt
  TABS="$TABS,$R/fam$f/refined_ptcl3D_print.txt"; SUBS="$SUBS,$R/expand/family${f}_subset.txt"; done
echo "== $M ($ST): family refine3D_states poses"; $PY $W/leaf_pose_bias.py $W/sim$M/all.txt ${TABS#,} $ST ${SUBS#,} 2>&1 | grep -v Warning
echo "== $M: initial consensus poses"; $PY $W/leaf_pose_bias.py $W/sim$M/all.txt $W/lat_$M/ptcl3D_r0.txt $ST 2>&1 | grep -v Warning

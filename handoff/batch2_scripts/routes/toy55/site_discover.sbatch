#!/bin/bash
#SBATCH -J klem_site_disc
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=160G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# site discovery on one node (compositional-first step, steps 1-2): many small view-balanced groups from the node's flex latent
# (site_groups.py, K = N / SIZE), their maps (reconstruct3D at the node's poses), the groups' cross-half variance -> significant
# reproducible variability components (site_components.py: cross-half eigen-decomposition, null from the negative eigenvalues), optional truth regions. env PROJ (project with the node's poses), FLEX (flex_pca_coordinates.txt),
# ROWS (per-row node membership; default all rows), GT, TRUTH ("name:A.mrc:B.mrc ..."), SIZE (8000), MSKDIAM (220), TAG (output name)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem
PROJ=${PROJ:?}; FLEX=${FLEX:?}; MD=${MSKDIAM:-220}; O=$K_/site_${TAG:?}; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') site $TAG] $*"; }
[ -s $O/ptcl3D.txt ] || simple_exec prg=print_project_field projfile=$PROJ oritype=ptcl3D mkdir=no > $O/ptcl3D.txt 2> /dev/null
NR=$(grep -c "e1=" $O/ptcl3D.txt); [ -n "${ROWS:-}" ] && cp -f $ROWS $O/rows.txt || awk -v n=$NR 'BEGIN {for (i = 1; i <= n; i++) print 1}' > $O/rows.txt
GX=""; [ -n "${GT:-}" ] && GX="--gt $GT"
[ -s $O/groups.txt ] || $PY $T/site_groups.py --oris $O/ptcl3D.txt --rows $O/rows.txt --flex $FLEX --out $O --size ${SIZE:-8000} $GX 2>&1 | grep -v Warning
KG=$(awk '$1 > 0' $O/groups.txt | sort -n | tail -1); log "$KG groups"
if [ -z "$(ls $O/rec/*_reconstruct3D/recvol_state$(printf %02d $KG)_even_unfil.mrc 2>/dev/null)" ]; then mkdir -p $O/rec && cd $O/rec && cp -f $PROJ g.simple
  simple_exec prg=selection projfile=g.simple oritype=ptcl3D infile=$O/groups.txt mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=g.simple qsys_name=local nparts=4 nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=g.simple pgrp=c1 mskdiam=$MD nstates=$KG nparts=4 nthr=8 > rec.log 2>&1 || { log "reconstruct FAILED: $(grep -a -m2 'ERROR\|THROW' rec.log | cut -c1-150)"; exit 1; }
  cd $O; log "groups reconstructed"; fi
RR=$(ls -td $O/rec/*_reconstruct3D | head -1); : > $O/list.txt
for s in $(seq 1 $KG); do b=$RR/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] && echo "g$s ${b}_even_unfil.mrc ${b}_odd_unfil.mrc $(awk -v s=$s '$1 == s' $O/groups.txt | wc -l)" >> $O/list.txt; done
$PY $T/site_components.py $O/comps --list $O/list.txt --mskdiam $MD ${TRUTH:+--truth $TRUTH} 2>&1 | grep -v Warning
log DONE

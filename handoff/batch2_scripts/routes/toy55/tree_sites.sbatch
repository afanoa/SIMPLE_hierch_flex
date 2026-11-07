#!/bin/bash
#SBATCH -J klem_tree_sites
#SBATCH -p csbdevel
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# site layer inside a finished tree (tree_sites.py), then scoring (TRUTH) or pole-map readouts (READOUT); env RT OUT H MRA MRB CP CTFABS MSKDIAM NODES NSN GATE PARS PARF TRUTH TMAPS TMASKS ORACLE EXTRALABELS MAPN
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-csbdevel} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:?}; R=$K_/rflex_$RT; O=${OUT:?}; MD=${MSKDIAM:-220}; mkdir -p $O; cd $O
$PY -u $T/tree_sites.py --tree $R --out $O --h ${H:?} --cache ${CP:?} ${MRA:+--mra $MRA} ${MRB:+--mrb $MRB} --ctfabs ${CTFABS:-1} --mskdiam $MD ${NODES:+--nodes $NODES} --nsn ${NSN:-3} --gate ${GATE:-0.8} \
  --pars ${PARS:-4} --parf ${PARF:-2} ${TMASKS:+--tmasks $TMASKS} || { echo "tree_sites FAILED"; exit 1; }
if [ -n "${TRUTH:-}" ] && [ -s $O/sites_final.txt ]; then
  MX=""; [ -n "${TMAPS:-}" ] && MX="--maps $TMAPS --masks $TMASKS"
  $PY $T/score_tree_sites.py --tree $R --truth $TRUTH $MX --minn 500 --json $O/tree.json > $O/tree.txt 2>&1
  $PY $T/score_tree_sites.py --tree $R --truth $TRUTH --labels $O/hybrid_labels.txt $MX --minn 500 --json $O/treesites.json > $O/treesites.txt 2>&1
  $PY $T/compare_tree_sites.py tree:$O/tree.json treesites:$O/treesites.json ${ORACLE:+--oracle $ORACLE} 2>&1 | grep -v Warning
  FL=$(awk '!/^#/ {printf "%s=%s/fits/%s@%s/flex ", $1, o, $1, $2}' o=$O $O/sites_final.txt); EM=$(awk '!/^#/ {printf "%s_lo=%s/fits/%s@%s/lo %s_hi=%s/fits/%s@%s/hi ", $1, o, $1, $3, $1, o, $1, $3}' o=$O $O/sites_final.txt)
  $PY $T/score_site_flex.py --truth $TRUTH --rows $R/rows_r.txt --flex $FL --em $EM --labels tree=$R/final_labels.txt treesites=$O/hybrid_labels.txt ${EXTRALABELS:-} --json $O/score_flex.json 2>&1 | grep -v Warning | tee $O/score_flex.txt; fi
if [ "${READOUT:-0}" = 1 ] && [ -s $O/sites_final.txt ]; then
  CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc; V=$K_/view_rflex/treesites_$(basename $O); : > $O/spec_top.txt
  rec(){ local d=$1 lab=$2 src=$3; mkdir -p $d; [ -n "$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ] && return 0; cd $d && cp -f $src r.simple
    simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$lab mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=slurm nparts=8 nthr=8 mkdir=no > up.log 2>&1 \
      && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=1 nparts=8 nthr=8 > rec.log 2>&1 || echo "reconstruct $d FAILED"; cd $O; }
  while read nm home lev rest; do for pole in lo hi; do for N in ${MAPN:-2500 10000}; do f=$O/fits/$nm@$lev/$pole
    [ -s $f/top$N.txt ] || $PY -c "import numpy as np; p = np.loadtxt('$f/post_master.txt'); l = np.zeros(len(p), int); l[np.argsort(-p)[:$N]] = 1; l[p < 0] = 0; np.savetxt('$f/top$N.txt', l, fmt='%d')"
    rec $O/rec_${nm}_${pole}_top$N $f/top$N.txt $R/node_$lev/n.simple; RD=$(ls -td $O/rec_${nm}_${pole}_top$N/*_reconstruct3D 2>/dev/null | head -1)
    [ -n "$RD" ] && echo "${nm}_${pole}_top$N $RD/recvol_state01_even_unfil.mrc $RD/recvol_state01_odd_unfil.mrc" >> $O/spec_top.txt; done; done; done < <(grep -v '^#' $O/sites_final.txt)
  $PY $T/make_view.py $V $O/spec_top.txt --box 160 > /dev/null 2>&1; $PY $T/arm_readout.py $CONS $CL $O/arm_top.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec_top.txt) 2>&1 | grep -v Warning | tail -1
  grep -v '^#' $O/arm_top.txt | awk '{n = $1; sub(/.*\//, "", n); printf "%s closed-arm %s open-arm %s %s\n", n, $2, $3, $4}'; echo "maps (lp 10, unsharpened, full maps only): $V"; fi
echo "[$(date '+%F %T') tree sites $RT] DONE"

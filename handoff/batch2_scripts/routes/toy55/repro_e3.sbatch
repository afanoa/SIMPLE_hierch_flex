#!/bin/bash
#SBATCH -J klem_repro_e3
#SBATCH -p csbdevel
#SBATCH -c 4
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# deterministic seeding test for node f692st r1: seeds = 8 quantile bins along the first principal axis of one run's focused flex latent
# (LATENT = that run's flex_pca_coordinates.txt; nothing from other runs), references reconstructed from the seeds, then rflex13's stage A
# (same node project n.simple, same mask rsmask.mrc, same keys, focus3). Two jobs (script latent, port latent) differ only in the seeds.
# Writes rnasep_klem/repro_r1/e3_<TAG>/{seeds,stageA_states}.txt. env TAG, LATENT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=4
PY=/mnt/nasapps/production/python/3.12/bin/python3; K_=/mnt/beegfs/elmlund/afan/rnasep_klem; D=$K_/rflex_f692st/node_r1; O=$K_/repro_r1/e3_${TAG:?}; MD=232
mkdir -p $O/refs $O/stageA; log(){ echo "[$(date '+%F %T') repro_e3 $TAG] $*"; }
st(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s = 0; for (i = 1; i <= NF; i++) if ($i ~ /^state=/) {split($i, a, "="); s = a[2] + 0}; print s}'; }
[ -s $O/seeds.txt ] || $PY - <<PYEOF || { log "seeds FAILED"; exit 1; }
import numpy as np
F = np.loadtxt('${LATENT:?}', comments='#', ndmin=2); p = F[:, 0].astype(int); Z = F[:, 5:] - F[:, 5:].mean(0)
w, V = np.linalg.eigh(np.cov(Z.T)); pc1 = Z @ V[:, -1]; q = np.quantile(pc1, np.linspace(0, 1, 9)[1:-1]); b = 1 + np.searchsorted(q, pc1)
lab = np.zeros(sum(1 for _ in open('$D/refined_f.txt')), int); lab[p - 1] = b
np.savetxt('$O/seeds.txt', lab, fmt='%d'); print('PC1 share %.2f; seed sizes %s' % (w[-1] / w.sum(), np.bincount(b)[1:]))
PYEOF
if [ -z "$(ls $O/refs/*_reconstruct3D/recvol_state08.mrc 2>/dev/null)" ]; then cd $O/refs && cp -f $D/n.simple s.simple
  simple_exec prg=selection projfile=s.simple oritype=ptcl3D infile=$O/seeds.txt mkdir=no > sel.log 2>&1 && \
  simple_exec prg=reconstruct3D projfile=s.simple pgrp=c1 mskdiam=$MD nstates=8 nparts=8 nthr=8 > rec.log 2>&1 || { log "seed references FAILED"; exit 1; }; fi
RR=$(ls -td $O/refs/*_reconstruct3D | head -1); FV=""; for j in $(seq 1 8); do FV="$FV vol$j=$RR/recvol_state$(printf %02d $j).mrc"; done
if [ -z "$(ls $O/stageA/*_refine3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ]; then cd $O/stageA && cp -f $D/n.simple r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$O/seeds.txt mkdir=no > rsel.log 2>&1 || { log "selection FAILED"; exit 1; }
  log "stage A refine3D from the deterministic seeds"
  simple_exec prg=refine3D projfile=r.simple $FV nstates=8 refine=prob_neigh prob_neigh_mode=geom pgrp=c1 mskdiam=$MD lp=15 maxits=3 center=no rec_backend=pcg \
    pcg_mskfile=$D/rsmask.mrc pcg_focus=yes nspace=2500 box_crop=64 smpd_crop=3.66000000 nparts=8 nthr=8 > ref.log 2>&1 || { log "stage A FAILED: $(grep -a -m2 ERROR ref.log | cut -c1-150)"; exit 1; }; fi
RDA=$(ls -td $O/stageA/*_refine3D | head -1); st $RDA/r.simple > $O/stageA_states.txt; log "stage A done: $(awk '$1 > 0' $O/stageA_states.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"
log DONE

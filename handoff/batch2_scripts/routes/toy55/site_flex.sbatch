#!/bin/bash
#SBATCH -J klem_site_flex
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# focused flex per site: flex_pca pcg_mskfile=<site mask> on ROWS of SRC for MASKS 'name=mask ..'; env RT NODE SRC ROWS MEANREC OUT K LPF NP PARL MSKDIAM TRUTH EMS LABELS REFS READOUT TOPN
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=${SNAP:-$HOME/SIMPLE_focus/snap/focus3}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=${QPART:-csbdevel} OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; RT=${RT:?}; R=$K_/rflex_$RT; ND=${NODE:-r}
SRC=${SRC:-$R/node_$ND/n.simple}; ROWS=${ROWS:-$R/rows_$ND.txt}; MD=${MSKDIAM:-220}; O=${OUT:?}; V=$K_/view_rflex/sflex_$(basename $O); mkdir -p $O; cd $O
log(){ echo "[$(date '+%F %T') site flex $RT] $*"; }; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1); CONS=$T/maps160/toy_consensus.mrc; NMR=$(grep -c . $ROWS)
rec(){ local d=$1 lab=$2; mkdir -p $d; [ -n "$(ls $d/*_reconstruct3D/recvol_state01_even_unfil.mrc 2>/dev/null)" ] && return 0; cd $d && cp -f $SRC r.simple
  simple_exec prg=selection projfile=r.simple oritype=ptcl3D infile=$lab mkdir=no > sel.log 2>&1 && simple_exec prg=update_project projfile=r.simple qsys_name=slurm nparts=${NP:-8} nthr=8 mkdir=no > up.log 2>&1 \
    && simple_exec prg=reconstruct3D projfile=r.simple pgrp=c1 mskdiam=$MD nstates=1 nparts=${NP:-8} nthr=8 > rec.log 2>&1 || log "reconstruct $d FAILED: $(grep -a -m2 ERROR rec.log | cut -c1-150)"; cd $O; }
if [ "${MEANREC:-0}" = 1 ]; then rec $O/mean $ROWS; MEAN=$(ls -td $O/mean/*_reconstruct3D | head -1)/recvol_state01.mrc; else MEAN=$(ls -td $R/node_$ND/mean/*_reconstruct3D | head -1)/recvol_state01.mrc; fi
[ -s $MEAN ] || { log "no mean map $MEAN"; exit 1; }
sflex(){ local nm=$1 mk=$2 D=$O/$1; mkdir -p $D; [ -n "$(ls $D/*_flex_pca/flex_pca_coordinates.txt 2>/dev/null)" ] && return 0; cd $D && cp -f $SRC f.simple
  simple_exec prg=selection projfile=f.simple oritype=ptcl3D infile=$ROWS mkdir=no > fsel.log 2>&1 || { log "site $nm selection FAILED"; return 1; }
  simple_exec prg=flex_pca projfile=f.simple vol1=$MEAN mskdiam=$MD lp=${LPF:-15} npreimages=${K:-8} rec_states=yes min_state_frac=0 umap=no n_probe_iters=${NPROBE:-2} pcg_mskfile=$mk nparts=${NP:-8} nthr=8 > flex.log 2>&1 \
    || { log "site $nm flex_pca FAILED: $(grep -a -m2 'ERROR\|THROW' flex.log | cut -c1-200)"; return 1; }
  log "site $nm: focused flex done ($(grep -a 'Execution time' flex.log | tail -1 | awk '{print $(NF-1)}') s), state sizes $(awk '!/^#/ {print $3}' $(ls -td $D/*_flex_pca | head -1)/flex_pca_coordinates.txt | sort -n | uniq -c | awk '{printf "%d ", $1}')"; }
log "$(echo ${MASKS:?} | wc -w) sites: focused flex on $(awk '$1 > 0' $ROWS | wc -l) rows of $SRC (mean $MEAN), ${PARL:-1} at a time"
for nm_mk in $MASKS; do ( sflex ${nm_mk%%=*} ${nm_mk#*=} ) & while [ $(jobs -rp | wc -l) -ge ${PARL:-1} ]; do sleep 30; done; done; wait
FL=""; for nm_mk in $MASKS; do FL="$FL ${nm_mk%%=*}=$O/${nm_mk%%=*}"; done
if [ -n "${TRUTH:-}" ]; then $PY $T/score_site_flex.py --truth $TRUTH --rows $ROWS --flex $FL ${REFS:+--ref $REFS} ${EMS:+--em $EMS} ${LABELS:+--labels $LABELS} --json $O/score.json 2>&1 | grep -v Warning | tee $O/score.txt; fi
if [ "${READOUT:-0}" = 1 ]; then
  : > $O/spec_states.txt; for e in $FL ${REFS:-}; do nm=${e%%=*}; FD=$(ls -td ${e#*=}/*_flex_pca | head -1)
    for s in $(seq 1 ${K:-8}); do b=$FD/recvol_state$(printf %02d $s); [ -s ${b}_even_unfil.mrc ] && echo "${nm}_state$s ${b}_even_unfil.mrc ${b}_odd_unfil.mrc" >> $O/spec_states.txt; done; done
  $PY $T/make_view.py $V $O/spec_states.txt --box 160 > /dev/null 2>&1; $PY $T/arm_readout.py $CONS $CL $O/arm_states.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec_states.txt) 2>&1 | grep -v Warning | tail -1
  $PY - <<PYEOF || { log "top-N selection FAILED"; exit 1; }
import numpy as np, glob, os
arm = {os.path.basename(w[0])[:-4]: (float(w[1]), float(w[2])) for w in (l.split() for l in open('$O/arm_states.txt')) if w and not w[0].startswith('#')}
for e in '$FL ${REFS:-}'.split():
    nm, d = e.split('=', 1); fd = sorted(glob.glob(d + '/*_flex_pca'), key=os.path.getmtime)[-1]; C = np.loadtxt(fd + '/flex_pca_coordinates.txt', comments='#')
    r, lab, Z = C[:, 0].astype(int) - 1, C[:, 2].astype(int), C[:, 5:]; st = [k for k in arm if k.startswith(nm + '_state')]; s = max(st, key=lambda k: arm[k][0]); sid = int(s.split('_state')[1])
    u = Z[lab == sid].mean(0) - Z.mean(0); o = np.argsort(-((Z - Z.mean(0)) @ u)); os.makedirs('$O/top_' + nm, exist_ok=True)
    for N in [int(x) for x in '${TOPN:-2500,5000,10000}'.split(',')]:
        l = np.zeros($NMR, int); l[r[o[:N]]] = 1; np.savetxt('$O/top_%s/top%d.txt' % (nm, N), l, fmt='%d')
    print('%s: states (closed / open, size) %s; most closed state %d -> top-N along its direction' % (nm, ' '.join('%d:%.2f/%.2f,%d' % (int(k.split('_state')[1]), arm[k][0], arm[k][1], (lab == int(k.split('_state')[1])).sum()) for k in sorted(st, key=lambda k: int(k.split('_state')[1]))), sid))
PYEOF
  : > $O/spec_top.txt; for e in $FL ${REFS:-}; do nm=${e%%=*}; for N in $(echo ${TOPN:-2500,5000,10000} | tr ',' ' '); do rec $O/rec_${nm}_top$N $O/top_$nm/top$N.txt
    RD=$(ls -td $O/rec_${nm}_top$N/*_reconstruct3D 2>/dev/null | head -1); [ -n "$RD" ] && echo "${nm}_top$N $RD/recvol_state01_even_unfil.mrc $RD/recvol_state01_odd_unfil.mrc" >> $O/spec_top.txt; done; done
  $PY $T/make_view.py $V $O/spec_top.txt --box 160 > /dev/null 2>&1; $PY $T/arm_readout.py $CONS $CL $O/arm_top.txt $(awk -v v=$V '{print v "/" $1 ".mrc"}' $O/spec_top.txt) 2>&1 | grep -v Warning | tail -1
  grep -v '^#' $O/arm_top.txt | awk '{n = $1; sub(/.*\//, "", n); printf "%s closed-arm %s open-arm %s %s\n", n, $2, $3, $4}'
  $PY $T/radial_extent.py $CONS $CL --maps $(grep -v '^#' $O/arm_top.txt | sort -k2,2gr | head -6 | awk '{n = $1; sub(/.*\//, "", n); sub(/\.mrc$/, "", n); printf "%s:%s ", n, $1}') 2>&1 | grep -v Warning
  log "maps (lp 10, unsharpened, full maps only): $V"; fi
log DONE

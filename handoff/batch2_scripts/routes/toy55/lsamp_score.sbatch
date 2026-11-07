#!/bin/bash
#SBATCH -J klem_lsamp_score
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# landscape sampling, step 1: score space of a finished rflex tree (lsamp_score.py): every tree particle vs every leaf map, cross-half, at its
# own leaf-node pose. Builds rflex_<RT>/lsamp/poses_combined.txt (each master row's line from its leaf's node project) unless POSES is given.
# env RT (rflex_<RT>, read-only) or TREE (any dir with leaves_list.txt + final_labels.txt), OUT (rnasep_klem/lsamp/<name>), SIM + MR (half stacks S_A,S_B + meta_rows MRA,MRB, or one stack + MR=full),
# CP (bench cache prefix), CTFABS (0 synthetic / 1 real phase-flipped), EMLIM (20), POSES (optional master-order pose table)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done
B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; R=${TREE:-$K_/rflex_${RT:?}}; RT=${RT:-$(basename $R)}; O=${OUT:?}
mkdir -p $O; log(){ echo "[$(date '+%F %T') lsamp_score $RT] $*"; }
P=${POSES:-$O/poses_combined.txt}
if [ ! -s $P ]; then
  for pj in $(awk '{print $2}' $R/leaves_list.txt | xargs -n1 dirname | sort -u); do f=$(ls $pj/*.simple | head -1); k=$(echo $pj | md5sum | cut -c1-12)
    [ -s $O/dump_$k.txt ] || simple_exec prg=print_project_field projfile=$f oritype=ptcl3D mkdir=no > $O/dump_$k.txt 2> /dev/null; log "dumped $f ($(grep -c . $O/dump_$k.txt) rows)"; done
  $PY - <<PYEOF
import hashlib, os, numpy as np
L = [l.split() for l in open('$R/leaves_list.txt') if l.strip()]; lab = np.loadtxt('$R/final_labels.txt', dtype=int); D = {}
for l in L:
    k = hashlib.md5((os.path.dirname(l[1]) + '\n').encode()).hexdigest()[:12]
    if k not in D: D[k] = open('$O/dump_%s.txt' % k).read().splitlines()
fb = D[next(iter(D))]; out = []
for r, s in enumerate(lab):   # tree rows: the leaf's node pose; others: any node's line (not scored)
    if s > 0: k = hashlib.md5((os.path.dirname(L[s - 1][1]) + '\n').encode()).hexdigest()[:12]; out.append(D[k][r])
    else: out.append(fb[r])
open('$P', 'w').write('\n'.join(out) + '\n'); print('pose table: %d rows (%d in the tree)' % (len(out), (lab > 0).sum()))
PYEOF
fi
log "scoring: leaves $R/leaves_list.txt, labels $R/final_labels.txt, poses $P, em_lim ${EMLIM:-20}"
$PY -u $T/lsamp_score.py --sim ${SIM:?} --mr ${MR:?} --cache ${CP:?} --poses $P --leaves $R/leaves_list.txt --labels $R/final_labels.txt --out $O \
    --em_lim ${EMLIM:-20} --ctfabs ${CTFABS:-0} --mskdiam ${MSKDIAM:-230} 2>&1 | grep -v "Warning\|warn("
log "DONE"

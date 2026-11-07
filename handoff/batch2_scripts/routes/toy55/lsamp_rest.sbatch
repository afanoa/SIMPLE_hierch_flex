#!/bin/bash
#SBATCH -J klem_lsamp_rest
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 01:30:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# how the rest state enters the pole-pair landscapes: per ITEM (NAME:NODEDIR:ROLES) the score space against the node's stage-B states
# (rest included) at its stage-B poses, the pair coordinates (lsamp_pair_coords.py; optional particle subsets by state), and equal-N held-out
# bins per pair tag (lsamp_repro_bins3.py; each coordinate has pole a low). env ITEMS (space-separated), PAIRS (TAG:A:B[:S+S],..), OUT
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64 OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; T=$A/routes/toy55; K_=$A/rnasep_klem; O=${OUT:?}; mkdir -p $O
H=$T/simrnasep692; MR=${H}_A/meta_rows.txt,${H}_B/meta_rows.txt; CP=$K_/lsamp/caches/rnasep692; CL=$(ls $T/maps160/docked/conf18_*_m160.mrc | head -1)
F=$K_/factorial/f692st_r1_polesets; REF=$F/script_z1lo/post_master.txt,$F/script_z1hi/post_master.txt,$F/script_z3hi/post_master.txt
for it in ${ITEMS:?}; do IFS=: read -r nm N roles <<< "$it"; pairs=${PAIRS:?}; D=$O/${nm}_tree; mkdir -p $D
  if [ ! -s $D/leaves_list.txt ]; then
    ( for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$"); do unset $v; done; B=$HOME/SIMPLE_focus/snap/focus3; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local
      cd $D && simple_exec prg=print_project_field projfile=$N/stageB/b.simple oritype=ptcl3D mkdir=no > dump.txt 2> /dev/null )
    $PY - <<PYEOF || { echo "$nm: tree setup FAILED"; continue; }
import numpy as np, glob
def kv(l): return dict(t.split('=', 1) for t in l.split() if '=' in t)
L = open('$D/dump.txt').read().splitlines(); r = np.loadtxt('$N/rows.txt', dtype=int) - 1; st = np.array([int(float(kv(L[i]).get('state', 0))) for i in r])
lab = np.zeros(len(L), int); lab[r] = st; np.savetxt('$D/final_labels.txt', lab, fmt='%d'); K = len(glob.glob('$N/stageB/recvol_state*_even_unfil.mrc'))
cnt = np.bincount(st, minlength=K + 1)[1:]; print('$nm: %d node particles, %d stage-B states %s' % (len(r), K, ' '.join(map(str, cnt))))
open('$D/leaves_list.txt', 'w').write(''.join('n_s%d $N/stageB/recvol_state%02d_even_unfil.mrc $N/stageB/recvol_state%02d_odd_unfil.mrc %d\n' % (k, k, k, cnt[k - 1]) for k in range(1, K + 1)))
PYEOF
  fi
  [ -s $O/${nm}_ps/meta.npz ] || { TREE=$D OUT=$O/${nm}_ps SIM=${H}_A,${H}_B MR=$MR CP=$CP CTFABS=1 EMLIM=12 MSKDIAM=232 POSES=$D/dump.txt bash $T/lsamp_score2.sbatch 2>&1 | grep "own leaf ranks\|DONE\|Error\|Traceback"; }
  RA=""; [ "$roles" = auto ] && RA="--roles auto --ref_post $REF"
  $PY -u $T/lsamp_pair_coords.py --scores $O/${nm}_ps --name $nm --out $O --pairs $pairs $RA 2>&1 | grep -v "Warning\|warn("
  for tg in $(echo $pairs | tr ',' '\n' | cut -d: -f1); do
    $PY -u $T/lsamp_repro_bins3.py --name $nm --rundir $O/$tg/$nm --scores $O/${nm}_ps --poses $D/dump.txt --ref_rundir $O/$tg/$nm --sim ${H}_A,${H}_B --mr $MR --cache $CP --ctfabs 1 --mskdiam 232 --out $O/$tg/bins 2>&1 | grep -v "Warning\|warn(\|map data"
    for f in A B; do $PY $T/arm_readout.py $T/maps160/toy_consensus.mrc $CL $O/$tg/bins/$nm/arm_bins_$f.txt $O/$tg/bins/$nm/bins_$f 2>&1 | grep -v Warning > /dev/null; done; done
done
echo "[$(date '+%F %T')] DONE"

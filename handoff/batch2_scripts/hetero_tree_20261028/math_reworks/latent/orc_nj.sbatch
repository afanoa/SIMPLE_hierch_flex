#!/bin/bash
#SBATCH -J orcnj
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%j.out
# Addendum to orc_run(2).sbatch on an existing oracle/<KEY>: the truth classes WITHOUT a junk reference (junk excluded from the
# maps, state 0 in the forced runs; assigned like everyone in the arms, starting label 1), at stage-B poses (F_njB, arm NB) and at
# true poses (F_njT, arm NT); then score_nj.txt with the arms already there. Env: KEY NODE.
set -uo pipefail
: ${KEY:?} ${NODE:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
SE=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$SE/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
export SIMPLE_PATH=$SE SIMPLE_SEED=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1
W=$L/oracle/$KEY
test -s $W/prep/lab_truth.txt || { echo "no $W/prep"; exit 4; }
cd $W
dump(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | grep 'e1=' > $2; }
BASE=$(grep -o 'simple_exec prg=refine3D .*' $NODE/stageB/hetero_step.sh | sed 's/ > hetero_step.log.*//' | sed -E "s/ vol[0-9]+=[^ ]+//g; s/ refine=[^ ]+/ refine=state_eval/; s/ prob_neigh_mode=[^ ]+//; s/ prob_beta=[^ ]+//; s/ nstates=[^ ]+//; s/ maxits=[^ ]+//; s/ projfile=[^ ]+//")
RBASE=$(grep -o 'simple_exec prg=reconstruct3D .*' $NODE/recB/hetero_step.sh | sed 's/ > hetero_step.log.*//' | sed -E "s/ nstates=[^ ]+//; s/ projfile=[^ ]+//")
# the label sets: nj_maps (junk 0), nj (junk 1), nj_map.txt (label class n), the forced prior of nj_maps
$PY - prep <<'EOF' || exit 2
import sys, os, numpy as np
p = sys.argv[1]
tm = np.loadtxt(os.path.join(p, 'truth_map.txt'), dtype=int, ndmin=2); lab = np.loadtxt(os.path.join(p, 'lab_truth.txt'), dtype=int)
jl = tm[tm[:, 1] == 0, 0]; keep = tm[tm[:, 1] != 0]
new = np.zeros(lab.max() + 1, int); new[keep[:, 0]] = np.arange(1, len(keep) + 1)
m = new[lab]; m[np.isin(lab, jl)] = 0; a = m.copy(); a[np.isin(lab, jl)] = 1
np.savetxt(os.path.join(p, 'lab_nj_maps.txt'), m, fmt='%d'); np.savetxt(os.path.join(p, 'lab_nj.txt'), a, fmt='%d')
np.savetxt(os.path.join(p, 'nj_map.txt'), np.c_[np.arange(1, len(keep) + 1), keep[:, 1], keep[:, 2]], fmt='%d')
K = len(keep)
with open(os.path.join(p, 'prior_nj_maps.txt'), 'w') as f:
    for i in np.nonzero(m > 0)[0]:
        f.write('%d %s\n' % (i + 1, ' '.join('0' if k == m[i] else '-1e4' for k in range(1, K + 1))))
print('no-junk label set: K %d, classes %s, junk rows %d' % (K, keep[:, 1].tolist(), (m == 0).sum()))
EOF
KN=$(sort -n prep/lab_nj.txt | tail -1)
mk(){ local out=$1 st=$2 po=${3:-}
  [ -s ${out%.simple}_dump.txt ] && return 0
  $PY $L/orc_patch.py b.simple $out check=prep/b_dump.txt states=$st ${po:+poses=$po} || exit 3
  dump $out ${out%.simple}_dump.txt
  $PY $L/orc_verify.py ${out%.simple}_dump.txt states=$st ${po:+poses=$po} || exit 3; }
mk pBnm.simple prep/lab_nj_maps.txt
mk pTnm.simple prep/lab_nj_maps.txt prep/poses_true.txt
mk pBn.simple prep/lab_nj.txt
mk pTn.simple prep/lab_nj.txt prep/poses_true.txt
seval(){ local D=$W/$1 PRJ=$W/$2 K=$3 VD=$4 MAXITS=$5 PRIOR=${6:-}
  [ -s $D/after.txt ] && { echo "$1: done before"; return 0; }
  mkdir -p $D/stderrout; cd $D; cp $PRJ s.simple
  local VOLS=$(for s in $(seq 1 $K); do printf 'vol%d=%s/recvol_state%02d.mrc ' $s $VD $s; done)
  local CMD="$SE/bin/$BASE projfile=s.simple nstates=$K maxits=$MAXITS $VOLS"
  printf '%s\nprior=%s\n' "$CMD" "${PRIOR:-none}" > seval_cmd.txt
  if [ -n "$PRIOR" ]; then export SIMPLE_STATE_PRIOR=$W/$PRIOR; else unset SIMPLE_STATE_PRIOR; fi
  /usr/bin/time -f "wall %e s" $CMD > seval.log 2> seval_time.txt; local rc=$?
  unset SIMPLE_STATE_PRIOR
  echo "$1: $2 K $K maps $VD maxits $MAXITS prior ${PRIOR:-none}: exit $rc, $(tail -1 seval_time.txt)"
  if [ $rc != 0 ]; then grep -a -m3 ERROR seval.log; cd $W; return $rc; fi
  dump s.simple after.txt; cd $W; }
recon(){ local D=$W/$1 PRJ=$W/$2 K=$3
  [ -s $D/recvol_state01.mrc ] && { echo "$1: done before"; return 0; }
  mkdir -p $D/stderrout; cd $D; cp $PRJ s.simple
  local CMD="$SE/bin/$RBASE projfile=s.simple nstates=$K"; echo "$CMD" > rec_cmd.txt
  /usr/bin/time -f "wall %e s" $CMD > rec.log 2> rec_time.txt; local rc=$?
  echo "$1: reconstruct3D $2 K $K: exit $rc, $(tail -1 rec_time.txt)"; cd $W; return $rc; }
kept(){ $PY - $W/$1/after.txt $W/prep/$2 <<'EOF'
import sys, numpy as np
s = np.array([int(float(dict(t.split('=', 1) for t in l.split() if '=' in t)['state'])) for l in open(sys.argv[1]) if 'e1=' in l])
ref = np.loadtxt(sys.argv[2], dtype=int)
print('%s: forced labels kept on %.5f of the rows' % (sys.argv[1], (s == ref).mean()))
EOF
}
recon R_njB pBnm.simple $KN || exit 5; seval F_njB pBnm.simple $KN $W/R_njB 2 prep/prior_nj_maps.txt || exit 5; kept F_njB lab_nj_maps.txt
recon R_njT pTnm.simple $KN || exit 5; seval F_njT pTnm.simple $KN $W/R_njT 2 prep/prior_nj_maps.txt || exit 5; kept F_njT lab_nj_maps.txt
seval A_NB pBn.simple $KN $W/F_njB 1 || exit 6
seval A_NT pTn.simple $KN $W/F_njT 1 || exit 6
ARMS="P:$W/A_P:maj TB:$W/A_TB:maj NB:$W/A_NB:tmap=$W/prep/nj_map.txt NT:$W/A_NT:tmap=$W/prep/nj_map.txt OB:$W/A_OB:truth OT:$W/A_OT:truth"
[ -s A_CB/after.txt ] && ARMS="$ARMS CB:$W/A_CB:coarse"
[ -s A_CT/after.txt ] && ARMS="$ARMS CT:$W/A_CT:coarse"
$PY $L/orc_score.py prep score_nj.txt $ARMS
echo "$KEY nj done $(date -Is)"

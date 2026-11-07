#!/bin/bash
#SBATCH -J orc2
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH --mail-type=NONE
#SBATCH --output=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent/logs/%x_%j.out
# orc_run.sbatch, resumable (a step whose output exists is skipped), plus the fixed-K corner: F_cT = the truth merged into stage B's
# states at true poses, arm CT (T x F_cT). So both factorials are complete: poses {B,T} x maps {stage B, coarse oracle} at stage B's
# K, and poses {B,T} x maps {oracle at B, oracle at T} at the truth's K. Env: KEY NODE TRUTH ALLTXT [MUNODE] [ORC_MINN ORC_MINFRAC:
# a child node's rare classes join junk as one 'other' reference].
set -uo pipefail
: ${KEY:?} ${NODE:?} ${TRUTH:?} ${ALLTXT:?}
L=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/math_reworks/latent
SE=/mnt/beegfs/elmlund/afan/hetero_tree_20261028/state_eval/install
PY=/mnt/nasapps/production/python/3.12/bin/python3
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$SE/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
export SIMPLE_PATH=$SE SIMPLE_SEED=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1
W=$L/oracle/$KEY
mkdir -p $W/prep; cd $W
dump(){ simple_exec prg=print_project_field projfile=$1 oritype=ptcl3D mkdir=no 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' | grep 'e1=' > $2; }
BASE=$(grep -o 'simple_exec prg=refine3D .*' $NODE/stageB/hetero_step.sh | sed 's/ > hetero_step.log.*//' | sed -E "s/ vol[0-9]+=[^ ]+//g; s/ refine=[^ ]+/ refine=state_eval/; s/ prob_neigh_mode=[^ ]+//; s/ prob_beta=[^ ]+//; s/ nstates=[^ ]+//; s/ maxits=[^ ]+//; s/ projfile=[^ ]+//")
RBASE=$(grep -o 'simple_exec prg=reconstruct3D .*' $NODE/recB/hetero_step.sh | sed 's/ > hetero_step.log.*//' | sed -E "s/ nstates=[^ ]+//; s/ projfile=[^ ]+//")
KB=$(ls $NODE/stageB/recvol_state[0-9][0-9].mrc | wc -l)
echo "$KEY: node $NODE, stage-B states $KB; ORC_MINN=${ORC_MINN:-} ORC_MINFRAC=${ORC_MINFRAC:-}"; echo "refine3D: $BASE"; echo "reconstruct3D: $RBASE"

if [ ! -s prep/ptcl_info.txt ]; then
  cp $NODE/stageB/b.simple b.simple
  dump b.simple prep/b_dump.txt
  $PY $L/orc_prep.py prep $NODE prep/b_dump.txt $TRUTH $ALLTXT || exit 2
fi
KT=$(sort -n prep/lab_truth.txt | tail -1); KC=$(sort -n prep/lab_coarse.txt | tail -1)
echo "K: stage B $KB, truth $KT, coarse $KC"

mk(){ local out=$1 st=$2 po=${3:-}
  [ -s ${out%.simple}_dump.txt ] && return 0
  $PY $L/orc_patch.py b.simple $out check=prep/b_dump.txt states=$st ${po:+poses=$po} || exit 3
  dump $out ${out%.simple}_dump.txt
  $PY $L/orc_verify.py ${out%.simple}_dump.txt states=$st ${po:+poses=$po} || exit 3; }
mk pBt.simple prep/lab_truth.txt
mk pBc.simple prep/lab_coarse.txt
mk pTs.simple prep/lab_stageB.txt prep/poses_true.txt
mk pTt.simple prep/lab_truth.txt prep/poses_true.txt
mk pTc.simple prep/lab_coarse.txt prep/poses_true.txt

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

MB=$NODE/stageB
seval F_ctl b.simple $KB $MB 2 prep/prior_stageB.txt || exit 5; kept F_ctl lab_stageB.txt
recon R_oB pBt.simple $KT || exit 5; seval F_oB pBt.simple $KT $W/R_oB 2 prep/prior_truth.txt || exit 5; kept F_oB lab_truth.txt
recon R_oT pTt.simple $KT || exit 5; seval F_oT pTt.simple $KT $W/R_oT 2 prep/prior_truth.txt || exit 5; kept F_oT lab_truth.txt
recon R_cB pBc.simple $KC || exit 5; seval F_cB pBc.simple $KC $W/R_cB 2 prep/prior_coarse.txt || exit 5; kept F_cB lab_coarse.txt
recon R_cT pTc.simple $KC || exit 5; seval F_cT pTc.simple $KC $W/R_cT 2 prep/prior_coarse.txt || exit 5; kept F_cT lab_coarse.txt

seval A_P   b.simple   $KB $MB      1 || exit 6
seval A_C   b.simple   $KB $W/F_ctl 1 || exit 6
seval A_OB  pBt.simple $KT $W/F_oB  1 || exit 6
seval A_CB  pBc.simple $KC $W/F_cB  1 || exit 6
seval A_TB  pTs.simple $KB $MB      1 || exit 6
seval A_OT  pTt.simple $KT $W/F_oT  1 || exit 6
seval A_OBT pBt.simple $KT $W/F_oT  1 || exit 6
seval A_OTB pTt.simple $KT $W/F_oB  1 || exit 6
seval A_CT  pTc.simple $KC $W/F_cT  1 || exit 6

$PY $L/orc_score.py prep score.txt P:$W/A_P:maj C:$W/A_C:maj CB:$W/A_CB:coarse TB:$W/A_TB:maj CT:$W/A_CT:coarse \
    OB:$W/A_OB:truth OBT:$W/A_OBT:truth OTB:$W/A_OTB:truth OT:$W/A_OT:truth
$PY $L/orc_wrong.py prep wrong.txt
if [ -n "${MUNODE:-}" ] && [ ! -s wrong_mu.txt ]; then
  mkdir -p mu; cp $MUNODE/stageB/b.simple mu/b.simple; dump mu/b.simple mu/b_dump.txt
  $PY $L/orc_prep.py mu $MUNODE mu/b_dump.txt $TRUTH $ALLTXT > /dev/null && \
  $PY $L/orc_wrong.py mu wrong_mu.txt $MUNODE/stageB/prob_state_evidence.txt
fi
echo "$KEY done $(date -Is)"

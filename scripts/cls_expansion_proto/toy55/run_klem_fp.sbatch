#!/bin/bash
#SBATCH -J toy10_klem_fp
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Kernel-latent EM (klem.py) on a toy: env TAG ("" or c), SIMDIR, PROJ, CSV, INIT (pca|random), FLEX (flex dir for the pca init), N1 N2 NFREE, LPS, SIGMAS, KTAG
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; export TOY_SMPD=1.464 OMP_NUM_THREADS=1
TAG=${TAG:-}; SIMDIR=${SIMDIR:-sim}; PJ=${PROJ:-proj}; OUT=$W/lat$TAG; T=$W/$SIMDIR; CSV=$W/${CSV:-rnasep_conformers10.csv}; STATES=$(awk -F, 'NR>1 {printf "%s%s", (NR>2?",":""), $2}' $CSV)
INIT=${INIT:-pca}; FLEX=${FLEX:-$OUT/flex12}; FD=$(ls -td $FLEX/*_flex_pca 2>/dev/null | head -1); KTAG=${KTAG:-}
if [ -z "${AXIS:-}" ]; then LOG=$(ls -t $W/toy10${TAG}_battery.*.out 2>/dev/null | head -1); AXIS=$(grep -o "truth axis (PC1 of the conformers, STATES order): .*" $LOG | head -1 | sed 's/.*: //'); fi
RP=$(ls -td $W/$PJ/*_refine3D | head -1); CONS=$RP/recvol_state01.mrc
D=$OUT/klem$KTAG; mkdir -p $D; echo "[$(date '+%F %T')] START klem$KTAG: init $INIT ($FD), out $D"
$PY $W/klem.py $T/all.mrcs ${TABLE:-$OUT/ptcl3D_r0.txt} $T/gt_particles.txt $D --cons $CONS --flexdir "$FD" --states $STATES --axis=$AXIS --n1 ${N1:-7} --n2 ${N2:-7} --n3 ${N3:-1} --nfree ${NFREE:-2} \
    --lps "${LPS:-20 20 20 14 14 14 10 10 10}" --sigmas "${SIGMAS:-1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7}" --nproc ${NPROC:-16} --init $INIT --nsub ${NSUB:-0} --nkseed ${NKSEED:-0} ${BAGPAR:+--bag_parents $BAGPAR} ${DECONV:+--deconv $DECONV} ${SPLIT:+--split_rounds $SPLIT} ${SUBSET:+--subset $SUBSET} ${FLEXPRIOR:+--flex_prior $FD:$FLEXPRIOR} --pnoise --balance 0.5 --devw 0.05 ${EXTRA:-} ${ADAPT:+--adapt $ADAPT} ${DEVW:+--devw $DEVW} ${FOCUS:+--focus $FOCUS} ${PMARG:+--pose_marg $PMARG} --pose_marg_n ${PMN:-7} ${SMARG:+--shift_marg $SMARG} ${PMM:+--pose_marg_mstep} ${TSTART:+--temp_start $TSTART} ${PNM:+--pnoise_m} ${TEMPF:+--temp_final $TEMPF} ${PNOISE:+--pnoise} ${PSCALE:+--pscale} ${MIXPROP:+--mixprop} ${BAL:+--balance $BAL} ${BALU:+--balance_until $BALU} ${WEIGHTS:+--weights $WEIGHTS} --mskdiam ${MSKD:-200} ${EXTRA:-} --kernel ${KERNEL:-lattice} --graph_after ${GAFTER:-3} --graph_nn ${GNN:-4} ${AVGS:+--avgs $AVGS --avg_mix ${AVGMIX:-0} --avg_mode ${AVGMODE:-unmix} --avg_prior ${AVGPRIOR:-1}} --nlayer ${NLAYER:-1} ${TIE:-} ${BAGS:+--bags $BAGS --bag_beta ${BAGBETA:-2}} --tag "klem$KTAG$TAG" 2>&1 | grep --line-buffered -v "^   [0-9]* / "
echo "== coverage $KTAG"; $PY $W/state_coverage.py $D/klem_latent.txt $T/gt_particles.txt --n1 ${N1:-7} --n2 ${N2:-7} --argmax --tag klem$KTAG 2>&1 | grep -v Warning | tail -16
echo "[$(date '+%F %T')] DONE"

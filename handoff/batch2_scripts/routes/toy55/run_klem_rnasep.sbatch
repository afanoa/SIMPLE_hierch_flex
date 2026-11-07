#!/bin/bash
#SBATCH -J rnasep_klem
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Kernel-latent EM on the real RNase P set (692k refined particles in 10 stacks, box 160 @ 1.464 A): env NSUB, INIT (random|npz), KTAG, N1 N2, LPS, SIGMAS, META (pose/CTF npz, default meta_all)
set -uo pipefail
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; export TOY_SMPD=1.464 OMP_NUM_THREADS=1
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
META=${META:-$A/rnasep_coupled/basis_span_all/meta_all.npz}; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
D=$A/rnasep_klem/klem${KTAG:-}; mkdir -p $D; echo "[$(date '+%F %T')] START rnasep klem${KTAG:-}: init ${INIT:-random}, nsub ${NSUB:-0}, out $D"
$PY $W/klem.py $STK none none $D --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 ${N1:-7} --n2 ${N2:-7} --nfree ${NFREE:-2} \
    --lps "${LPS:-20 20 20 14 14 14 10 10 10}" --sigmas "${SIGMAS:-1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7}" --nproc 16 --init ${INIT:-random} --nsub ${NSUB:-0} --ctf auto ${BAGPAR:+--bag_parents $BAGPAR} ${DECONV:+--deconv $DECONV} ${SPLIT:+--split_rounds $SPLIT} ${AVGS:+--avgs $AVGS --avg_mix ${AVGMIX:-0} --avg_mode ${AVGMODE:-unmix} --avg_prior ${AVGPRIOR:-1}} ${ADAPT:+--adapt $ADAPT} ${DEVW:+--devw $DEVW} ${FOCUS:+--focus $FOCUS} ${PMARG:+--pose_marg $PMARG} --pose_marg_n ${PMN:-7} ${SMARG:+--shift_marg $SMARG} ${PMM:+--pose_marg_mstep} ${TSTART:+--temp_start $TSTART} ${PNM:+--pnoise_m} ${TEMPF:+--temp_final $TEMPF} ${PNOISE:+--pnoise} ${PSCALE:+--pscale} ${MIXPROP:+--mixprop} ${BAL:+--balance $BAL} ${BALU:+--balance_until $BALU} --kernel ${KERNEL:-lattice} --graph_after ${GAFTER:-3} --graph_nn ${GNN:-4} ${BAGS:+--bags $BAGS --bag_beta ${BAGBETA:-2}} --tag "rnasep${KTAG:-}" 2>&1 | grep --line-buffered -v "^   [0-9]* / "
echo "[$(date '+%F %T')] DONE"

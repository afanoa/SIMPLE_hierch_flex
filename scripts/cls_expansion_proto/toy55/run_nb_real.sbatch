#!/bin/bash
#SBATCH -J rnasep_nb_real
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=96G
#SBATCH -t 06:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# Python parent board on a 100k random subset of the 692k apo RNase P (|CTF|, contrast scale, pnoise, balance) with
# the mixture-of-PPCA node basis (env NB: 0 = baseline, q > 0): is there a cross-half reproducible within-node residual
# mode in the real data? Readouts: node-basis reproducibility/amplitude per band, identification vs 10 and 55 maps.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; export OMP_NUM_THREADS=1 TOY_SMPD=1.464
META=${META:-$A/rnasep_coupled/basis_span_all/meta_all.npz}; CONS=$A/rnasep_states/25_refine3D_states/rec_final_state01.mrc; LAB=$A/rnasep_softem/states.txt
STK=$(ls /mnt/beegfs/elmlund/rnasep/22_denoise_project/denoise_raw_particles_part??.mrcs | sort | tr "\n" , | sed 's/,$//')
NB=${NB:-0}; TAG=${TAG:-nb$NB}; R=$A/rnasep_klem/nb_real; D=$R/board_$TAG; mkdir -p $R $D; log(){ echo "[$(date '+%F %T') nbreal_$TAG] $*"; }
SUBSET=${SUBSET:-$R/subset100k.txt}; LPS=${LPS:-20 20 20 14 14 14 10 10 10}; SIG=${SIG:-1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7}; FOCUS=${FOCUS:-}
if [ -n "${FAMSUB:-}" ] && [ ! -s $SUBSET ]; then $PY -c "import numpy as np; s=np.loadtxt('$FAMSUB',dtype=int); r=np.random.default_rng(7); np.savetxt('$SUBSET', np.sort(r.choice(s,min(100000,len(s)),replace=False)), fmt='%d')"; fi
[ -s $R/subset100k.txt ] || $PY -c "import numpy as np; s=np.loadtxt('$A/rnasep_klem/f692c/subset_all.txt',dtype=int); r=np.random.default_rng(7); np.savetxt('$R/subset100k.txt', np.sort(r.choice(s,100000,replace=False)), fmt='%d')"
log "python board on $(grep -c . $SUBSET) particles (meta $(basename $META), focus ${FOCUS:-none}), |CTF|, pscale, node basis q=$NB ${EXTRA:-}"; t0=$(date +%s)
$PY $W/klem.py $STK none none $D --meta $META --labels $LAB --cons $CONS --states 1,2,3,4 --axis=0,0,0,0 --n1 7 --n2 7 --nfree 2 --nproc 16 --init npz --nsub 0 --subset $SUBSET --ctf abs ${FOCUS:+--focus_mask $FOCUS} --pscale --pnoise --balance 0.5 --devw 0.05 \
    --lps "$LPS" --sigmas "$SIG" $([ $NB -gt 0 ] && echo "--node_basis $NB --node_basis_after 3") ${EXTRA:-} --tag "nbreal_$TAG" 2>&1 | grep --line-buffered -v "^   [0-9]* / " | grep "node basis\|final:\|states: cross-half\|CTF model\|contrast\|Error\|Traceback\|it [0-9]*:" | cut -c1-260
log "board done in $(( ($(date +%s) - t0) / 60 )) min"
echo "== identification of the cells vs the 10 deposited"; $PY $W/deposited_match.py $D/cells $W/maps160/docked $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -2
echo "== identification of the cells vs ALL 55 deposited apo"; $PY $W/deposited_match.py $D/cells $W/maps160/docked_all $W/maps160/toy_consensus.mrc 2>&1 | grep -v Warning | tail -2
log DONE

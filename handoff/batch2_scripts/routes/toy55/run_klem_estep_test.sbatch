#!/bin/bash
#SBATCH -J toy10_klem_estep_test
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 02:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Fortran klem E-step vs the Python reference dump on the toy: inputs (particle rows + halves, per-shell noise/scale)
# built from klem_estep_ref.npz, simple_exec prg=klem_estep run with the dump's calibrated conventions, with the
# fftshift phase origin on and off, and compared coefficient by coefficient (compare_estep.py).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=${BDIR:-$HOME/SIMPLE_klem2/build}; DONE=${DONEFILE:-make0.done}; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=local OMP_NUM_THREADS=16
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; D=$W/lat/klem_dump20k; T=$D/fortran_test; mkdir -p $T; log(){ echo "[$(date '+%F %T') ft] $*"; }
n=0; while [ $n -lt 90 ] && { [ ! -f $B/$DONE ] || [ ! -s $D/klem_estep_ref.npz ] || ! grep -q F200 <(/mnt/nasapps/production/python/3.12/bin/python3 -c "import numpy as np; print(' '.join(np.load('$D/klem_estep_ref.npz').keys()))" 2>/dev/null); }; do sleep 60; n=$((n+1)); done
cat $B/$DONE; grep -q "_EXIT=0" $B/$DONE || { log "build failed"; exit 1; }
cd $T; $PY - <<PYEOF
import numpy as np
Z = np.load('$D/klem_estep_ref.npz'); pinds = Z['pinds'].astype(int); eo = Z['eo'].astype(int); sig2 = Z['sig2']; SC = Z['SC']; sh = Z['shell'].astype(int); kmax = int(Z['kmax'])
np.savetxt('$T/pinds_eo.txt', np.c_[pinds, eo], fmt='%d')
with open('$T/sig2_scale.txt', 'w') as f:
    for s in range(kmax + 1):
        m = sh == s; f.write('%d %.8e %.8e\n' % (s, float(sig2[m][0]) if m.any() else 1., float(SC[m][0]) if m.any() else 1.))
print('conventions: ss %s fs %s conv %s sign %s smpd %s | lp %s kmax %d n3 %s K %s nband %d particles %d' % (Z['ss'], Z['fs'], Z['conv'], Z['sign'], Z['smpd'], Z['lp'], kmax, Z['n3'], Z['K'], len(sig2), len(pinds)))
open('$T/env.sh', 'w').write('export KLEM_SS=%g KLEM_FS=%g KLEM_CONV=%d KLEM_SIGN=%g\n' % (float(Z['ss']), float(Z['fs']), int(Z['conv']), float(Z['sign'])))
open('$T/args.txt', 'w').write('nstates=%d lp=%g\n' % (int(Z['K']), float(Z['lp'])))
PYEOF
source $T/env.sh; ARGS=$(cat $T/args.txt); PJ=$(ls -td $W/proj/*_refine3D | head -1)/proj.simple; export KLEM_PINDS=$T/pinds_eo.txt KLEM_NSAMP=200
for ph in 1 0; do export KLEM_PHORIG=$ph; log "simple_exec prg=klem_estep ($ARGS) phase origin $ph"
  simple_exec prg=klem_estep projfile=$PJ dir=$D/cells_half $ARGS infile=$T/sig2_scale.txt outfile=$T/estep_ph$ph.bin nthr=16 mkdir=no > $T/estep_ph$ph.log 2>&1 || { log "klem_estep FAILED (phorig $ph): $(grep -a -m3 'ERROR\|STOP\|Error' $T/estep_ph$ph.log | cut -c1-200)"; tail -5 $T/estep_ph$ph.log; continue; }
  grep -a "KLEM" $T/estep_ph$ph.log | cut -c1-160; echo "== comparison (phase origin $ph)"; $PY $W/compare_estep.py $D $T/estep_ph$ph.bin 2>&1 | grep -v Warning
done
log DONE

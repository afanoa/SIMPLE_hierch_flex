#!/bin/bash
#SBATCH -J rnasep_holo_maps3
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 08:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# Holo (RNase P RNA + rnpA, 1 mM Mg) conformer maps for the apo/holo x open/closed synthetic test: download 8 of the
# 21 deposited holo conformers (EMD-74272..74279), resample/mirror to our lattice (process_maps.py, as the apo maps),
# dock each into the apo consensus frame (dock_volpair vs the docked conformer 1), then a readout per holo map:
# docking cc, closed/open by the P12/P10.1 contact (correlation inside the hand mask with closed-18 vs the open mean,
# arm occupancy), and the protein excess (where the holo map exceeds the apo consensus outside the arm) written as a
# soft mask (rnpa_mask.mrc) = the truth region of the compositional difference.
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$PATH OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; cd $W; log(){ echo "[$(date '+%F %T') holo] $*"; }
CSV=$W/rnasep_holo8.csv; mkdir -p $W/emdb $W/maps160/holo_docked $W/maps160/holo
for id in $(tail -n +2 $CSV | cut -d, -f1 | sed 's/EMD-//'); do f=$W/emdb/emd_$id.map.gz
  if [ -s $f ] && gzip -t $f 2>/dev/null; then log "have $f (verified)"; continue; else rm -f $f; fi
  for try in 1 2 3 4; do wget -q -c -O $f https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-$id/map/emd_$id.map.gz && gzip -t $f 2>/dev/null && break; sleep 30; done
  gzip -t $f 2>/dev/null || { log "download of emd_$id.map.gz FAILED verification"; rm -f $f; }
  log "$(date +%T) emd_$id.map.gz $(stat -c %s $f 2>/dev/null)"
done
$PY $W/process_maps.py $CSV $W/emdb $W/maps160 || { log "process FAILED"; exit 1; }
REF=$A/rnasep_coupled/degenhardt/emd_74131_m160_docked.mrc
for m in $W/maps160/conf1??_m160.mrc; do b=$(basename $m .mrc)
  [ -s $W/maps160/holo_docked/$b.mrc ] && continue
  simple_exec prg=dock_volpair vol1=$REF vol2=$m smpd=1.464 lp=10 hp=120 mskdiam=200 trs=12 outvol=$W/maps160/holo_docked/$b.mrc outfile=$W/maps160/holo_docked/$b.txt nthr=8 mkdir=no > $W/maps160/holo_docked/$b.log 2>&1 || { log "dock FAILED $b"; continue; }
  echo "$b $(grep '^cc=' $W/maps160/holo_docked/$b.txt) $(grep '^euler=' $W/maps160/holo_docked/$b.txt) $(grep '^shift=' $W/maps160/holo_docked/$b.txt)"
done
$PY - <<PYEOF
import numpy as np, glob, sys, csv
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
from scipy.ndimage import gaussian_filter, binary_dilation
W = '$W'; SMPD = 1.464; n = 160; g = lambda v: gaussian_filter(v, 2.0)
C = {i: read_mrc(glob.glob(W + '/maps160/docked/conf%02d_*.mrc' % i)[0]) for i in (1, 10, 18, 27, 28, 37, 46, 48, 49, 51)}; cons = read_mrc(W + '/maps160/toy_consensus.mrc')
# the hand mask of the contact (as closed_focal_rnasep.sbatch) and the closed-arm occupancy definition (as run_klem_matched_rec)
Lc = g(C[18]); Lo = g(np.mean([C[i] for i in (27, 28, 48, 49)], 0)); diff = Lc - Lo; dens = Lc > Lc.mean() + 1.5 * Lc.std()
core = (diff > diff.mean() + 2.5 * diff.std()) & dens; m = binary_dilation(core, iterations=4); mk = gaussian_filter(m.astype(float), 2.0); mk = np.clip(mk / max(mk.max(), 1e-30), 0, 1)
closed = g(C[18] - cons); openm = g(np.mean([C[i] for i in (27, 28, 48, 49)], 0) - cons)
def fc(x, y): return float(np.sum(mk * x * y) / max(np.sqrt(np.sum(mk * x * x) * np.sum(mk * y * y)), 1e-30))
lp25 = lambda v: gaussian_filter(v, 2.5); L18 = lp25(C[18]); Lcons = lp25(cons); dens18 = L18 > L18.mean() + 2.0 * L18.std(); d18 = L18 - Lcons
arm = (d18 > d18.mean() + 3 * d18.std()) & dens18; core18 = dens18 & (np.abs(d18) < 0.5 * d18.std()); ref = L18[arm].mean() / L18[core18].mean()
def occ(v): Lv = lp25(v); return (Lv[arm].mean() / Lv[core18].mean()) / ref
print('apo references: focal score (corr closed - corr open inside the contact mask) and arm occupancy')
for i, v in C.items(): d = g(v - cons); print('  apo conformer %2d: focal %+.2f  arm occ %.2f' % (i, fc(d, closed) - fc(d, openm), occ(v)))
rows = list(csv.DictReader(open('$CSV'))); excess = []
print('holo maps docked into the apo frame:')
for r in rows:
    fs = glob.glob(W + '/maps160/holo_docked/conf%d_*.mrc' % int(r['conformer']))
    if not fs: print('  holo %s: not docked' % r['conformer']); continue
    v = read_mrc(fs[0]); cc = float(open(fs[0][:-4] + '.txt').read().split('cc=')[1].split()[0]); d = g(v - cons)
    # protein excess: positive deviation from the apo consensus outside the arm region, inside/near the density
    ex = np.where(d > 0, d, 0.); ex[binary_dilation(arm, iterations=6)] = 0.; excess.append(ex)
    print('  holo %s (%s, %s ptcls): dock cc %.3f | focal %+.2f | arm occ %.2f | positive excess outside the arm: %.1f%% of the consensus mass above 2 sd' % (r['conformer'], r['emdb_id'], r['images'], cc, fc(d, closed) - fc(d, openm), occ(v), 100 * ex[ex > 2 * d.std()].sum() / max(Lcons[Lcons > 0].sum(), 1e-30)))
if excess:
    E = np.mean(excess, 0); thr = E.mean() + 3 * E.std(); mp = E > thr; mp = binary_dilation(mp, iterations=2); soft = gaussian_filter(mp.astype(float), 2.0); soft = np.clip(soft / max(soft.max(), 1e-30), 0, 1)
    write_mrc(W + '/maps160/holo/rnpa_mask.mrc', soft.astype(np.float32), SMPD); zz, yy, xx = np.indices((n, n, n)); rr = np.hypot(np.hypot(xx - n / 2, yy - n / 2), zz - n / 2) * SMPD
    print('mean holo excess over the apo consensus (outside the arm): %d voxels above 3 sd, radius median %.0f A, max %.0f A -> maps160/holo/rnpa_mask.mrc' % ((E > thr).sum(), np.median(rr[E > thr]) if (E > thr).any() else 0, rr[E > thr].max() if (E > thr).any() else 0))
    write_mrc(W + '/maps160/holo/holo_mean_excess.mrc', E.astype(np.float32), SMPD)
PYEOF
log DONE

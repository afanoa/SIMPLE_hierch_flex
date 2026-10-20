#!/bin/bash
#SBATCH -J klem_core_align
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The paper's rigid-core alignment on f692e family 1: refine3D of the family's particles against the family parent with a
# CORE mask (envelope minus the moving regions), local search from the family-refined poses; then the 4 test states
# (39k / 15k / 5k / closed 4.9k) reconstructed at the core-aligned poses: native FSC, arm-region FSC, closed occupancy,
# versus the same states at the family-refined poses (consolidated_12_fullw_*).
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692e; O=$R/core_align; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') corealign] $*"; }
PS=$(ls -td $R/fam1/refine/*_refine3D_states | head -1)/p.simple; CL=$(ls $W/maps160/docked/conf18_*_m160.mrc | head -1)
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
P = read_mrc('$R/parents/fam1.mrc').astype(np.float32); n = P.shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2
lp = lambda v, res: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (1.464 / res) ** 2))).astype(np.float32)
env = lp((lp(P, 10.) > 0.1 * lp(P, 10.).max()).astype(np.float32), 12.); mov = lp(read_mrc('$R/fam1/regions_motion/all_moving.mrc').astype(np.float32), 12.); mov = np.clip(mov / max(mov.max(), 1e-12), 0, 1)
core = np.clip(env * (1 - mov), 0, 1).astype(np.float32); write_mrc('$O/core_mask.mrc', core, 1.464); write_mrc('$O/parent_core.mrc', (P * core).astype(np.float32), 1.464); print('core mask: envelope %.3f of the box, moving regions removed %.3f, core %.3f' % (env.mean(), (env * mov).mean(), core.mean()))
PYEOF
cp -f $PS p.simple
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s=0; for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); s=a[2]}; print (s>0)?1:0}' > $O/all1.txt   # single-state project: every family particle in state 1
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/all1.txt mkdir=no > sel0.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=40000 nparts=20 nthr=8 mkdir=no > up.log 2>&1
log "refine3D with pcg_mskfile = core mask (reference solved inside the core), local search, lp 6, maxits 6 (20 x 8)"
simple_exec prg=refine3D projfile=p.simple vol1=$R/parents/fam1.mrc pcg_mskfile=$O/core_mask.mrc refine=neigh lp=6 maxits=6 mskdiam=220 pgrp=c1 nparts=20 nthr=8   # PCG support mask = the core: every iteration's reference is solved inside the core only (rigid-core alignment) > refine.log 2>&1 || { log "refine3D FAILED: $(grep -a -m2 ERROR refine.log | cut -c1-200)"; exit 1; }
[ "$(grep -a -c "RESOLUTION @ FSC=0.143" $O/refine.log)" -ge 2 ] || { log "refine3D did not iterate: $(grep -a -m2 "ERROR" $O/refine.log | cut -c1-200)"; exit 1; }
RD=$(ls -td $O/*_refine3D | head -1); grep -a "RESOLUTION @ FSC=0.143" $O/refine.log | tail -1 | cut -c1-120; log "core-aligned project $RD/p.simple"
cd $O; cp -f $RD/p.simple states.simple
$PY - <<PYEOF
import numpy as np
st = np.loadtxt('$R/consolidated_12_fullw_region/states.txt', dtype=int); zs = np.loadtxt('$R/consolidated_12_fullw_zoom/states.txt', dtype=int)
f1 = [int(l.split()[0]) for l in open('$R/consolidated_12_fullw_region/groups.txt') if not l.startswith('#') and l.split()[2].startswith('f1')]; sizes = np.bincount(st)
pick = [max(f1, key=lambda s: sizes[s]), min(f1, key=lambda s: abs(sizes[s] - 15000)), min(f1, key=lambda s: abs(sizes[s] - 5000))]; out = np.zeros(len(st), int)
for r, s in enumerate(pick): out[st == s] = r + 1
out[zs == 1] = 4; np.savetxt('$O/states.txt', out, fmt='%d'); print('family-1 test states: ' + ' '.join('%d:%d' % (r + 1, (out == r + 1).sum()) for r in range(4)) + ' (leaves %s + closed zoom state)' % pick)
np.savetxt('$O/before_states.txt', pick, fmt='%d')
PYEOF
simple_exec prg=selection projfile=states.simple oritype=ptcl3D infile=states.txt mkdir=no > sel.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=reconstruct3D projfile=states.simple pgrp=c1 mskdiam=232 nstates=4 nparts=10 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED"; exit 1; }
RR=$(ls -td $O/*_reconstruct3D | head -1); rm -rf "${O:?}/maps"; mkdir -p $O/maps; cp -f $RR/recvol_state0?*.mrc $O/maps/; find "${RR:?}" -name "*_part*.mrc" -delete
$PY - <<PYEOF
import numpy as np, sys, os
sys.path.insert(0, '$W'); from model_class import read_mrc
def fsc_res(fe, fo, m=None, smpd=1.464):
    e, o = read_mrc(fe).astype(np.float32), read_mrc(fo).astype(np.float32); n = e.shape[0]; x = np.arange(n) - n // 2
    rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1) * (1 if m is None else m)
    E, O = np.fft.fftn(e * msk), np.fft.fftn(o * msk); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2)
    num = np.bincount(sh.ravel(), np.real(np.conj(E) * O).ravel()); d1 = np.bincount(sh.ravel(), np.abs(E).ravel() ** 2); d2 = np.bincount(sh.ravel(), np.abs(O).ravel() ** 2)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < 0.143)[0]; k = (below[0] + 1) if len(below) else n // 2 - 1; return n * smpd / k
n = 160; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2
lp = lambda v, res: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (1.464 / res) ** 2))).astype(np.float32)
arm = np.clip(lp(read_mrc('$R/fam1/regions_motion/region_1.mrc').astype(np.float32), 15.), 0, 1) ** 2
pk = [int(x) for x in open('$O/before_states.txt').read().split()]; before = {i + 1: '$R/consolidated_12_fullw_region/fam0/regionall/level1/maps/recvol_state%02d' % pk[i] for i in range(3)}; before[4] = '$R/consolidated_12_fullw_zoom/fam0/regionall/level1/maps/recvol_state01'
print('state  global FSC (family poses -> core-aligned)   arm-region FSC (family -> core-aligned)')
for s in range(1, 5):
    b, a = before[s], '$O/maps/recvol_state%02d' % s
    print('%d   %.1f -> %.1f   %.1f -> %.1f' % (s, fsc_res(b + '_even.mrc', b + '_odd.mrc'), fsc_res(a + '_even.mrc', a + '_odd.mrc'), fsc_res(b + '_even.mrc', b + '_odd.mrc', arm), fsc_res(a + '_even.mrc', a + '_odd.mrc', arm)))
PYEOF
echo "== closed-arm occupancy at the core-aligned poses"; $PY $W/arm_readout.py $W/maps160/toy_consensus.mrc $CL $O/arm.txt $O/maps 2>&1 | grep -v Warning | grep -v mask | grep "recvol_state0[1-4].mrc\|summary"
log DONE

#!/bin/bash
#SBATCH -J klem_core_align_fam
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# rigid-core alignment of family FAM of f692e (env FAM): refine3D of all its particles against the family parent with
# pcg_mskfile = the core (envelope minus the moving regions), local search from the family-refined poses; the aligned
# project is left at $R/core_align_f$FAM/aligned.simple for the leaf-set rebuild
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; R=$A/rnasep_klem/f692e; F=${FAM:?}; O=$R/core_align_f$F; mkdir -p $O; cd $O; log(){ echo "[$(date '+%F %T') corealign f$F] $*"; }
PS=$(ls -td $R/fam$F/refine/*_refine3D_states | head -1)/p.simple
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc, write_mrc
P = read_mrc('$R/parents/fam$F.mrc').astype(np.float32); n = P.shape[0]; f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); k2 = kx**2 + ky**2 + kz**2
lp = lambda v, res: np.real(np.fft.ifftn(np.fft.fftn(v) * np.exp(-3. * k2 / (1.464 / res) ** 2))).astype(np.float32)
env = lp((lp(P, 10.) > 0.1 * lp(P, 10.).max()).astype(np.float32), 12.); mov = lp(read_mrc('$R/fam$F/regions_motion/all_moving.mrc').astype(np.float32), 12.); mov = np.clip(mov / max(mov.max(), 1e-12), 0, 1)
core = np.clip(env * (1 - mov), 0, 1).astype(np.float32); write_mrc('$O/core_mask.mrc', core, 1.464); print('core mask: envelope %.3f, moving removed %.3f, core %.3f of the box' % (env.mean(), (env * mov).mean(), core.mean()))
PYEOF
cp -f $PS p.simple
simple_exec prg=print_project_field projfile=p.simple oritype=ptcl3D mkdir=no 2>/dev/null | awk '{s=0; for(i=1;i<=NF;i++) if($i ~ /^state=/){split($i,a,"="); s=a[2]}; print (s>0)?1:0}' > $O/all1.txt
simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$O/all1.txt mkdir=no > sel0.log 2>&1 || { log "selection FAILED"; exit 1; }
simple_exec prg=update_project projfile=p.simple qsys_name=slurm qsys_partition=csbdevel job_memory_per_task=32000 walltime=40000 nparts=20 nthr=8 mkdir=no > up.log 2>&1
log "refine3D with pcg_mskfile = core, refine=neigh lp 6 maxits 6, 20 x 8 on $(grep -c 1 $O/all1.txt) particles"
simple_exec prg=refine3D projfile=p.simple vol1=$R/parents/fam$F.mrc pcg_mskfile=$O/core_mask.mrc refine=neigh lp=6 maxits=6 mskdiam=220 pgrp=c1 nparts=20 nthr=8 > refine_f$F.log 2>&1
RD=$(ls -td $O/*_refine3D | head -1); NI=$(ls $RD/RESOLUTION_STATE01_ITER* 2>/dev/null | wc -l); [ "$NI" -ge 2 ] || { log "refine3D did not iterate: $(grep -a -m2 ERROR refine_f$F.log | cut -c1-200)"; exit 1; }
cp -f $RD/p.simple $O/aligned.simple; log "DONE: $NI iterations, aligned project $O/aligned.simple"

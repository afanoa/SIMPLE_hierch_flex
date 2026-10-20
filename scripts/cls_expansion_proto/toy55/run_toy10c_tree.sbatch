#!/bin/bash
#SBATCH -J toy10c_tree
#SBATCH -p csbdevel
#SBATCH -c 16
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/routes/toy55/%x.%j.out
# The current recipe on the COMPOSITIONAL toy (toy10c: the 10 conformers + 61 = deleted domain, 62 = added density,
# 12 states + junk): parent board (junk comps) -> two families -> children (focus, no free comps) -> k-means 10 ->
# refine3D_states per family (poses only) -> children again at the refined poses (feedback recipe, balanced) -> level 3
# (balanced, focus, 12->8 A). Readout at every stage: purity per state including 61 and 62. env: BAL (0.5), PM (3)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
BK=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$BK PATH=$BK/bin:$BK/scripts:$PATH TOY_SMPD=1.464
export TAG=c SIMDIR=simc CSV=rnasep_conformers10c.csv PROJ=projc PJ=projc
PY=/mnt/nasapps/production/python/3.12/bin/python3; W=/mnt/beegfs/elmlund/afan/routes/toy55; cd $W; L=$W/latc; S=$W/simc; GT=$S/gt_particles.txt; log(){ echo "[$(date '+%F %T') c] $*"; }
log "1. parent board (random init, junk comps)"; env KTAG=_junk INIT=random NSUB=0 NFREE=2 bash $W/run_klem.sbatch 2>&1 | grep "final:\|states: cross-half\|Error\|Traceback" | sed 's/^/[parent] /'
$PY $W/state_coverage.py $L/klem_junk/klem_latent.txt $GT --argmax --tag "toy10c parent" | tail -1
log "2. families and children at consensus poses (focus, no free comps)"; env PARENT=$L/klem_junk ETAG=_c NFAM=2 FOCUS=1 NFREE=0 bash $W/run_expand_toyc.sbatch 2>&1 | grep "family [0-9]:\|between\|final:\|MERGED\|best-state\|Error\|Traceback" | sed 's/^/[children] /'
log "3. per-family 10-state refinement (poses only)"; RT=$L/ptcl3D_realigned_c.txt
for j in 1 2; do
  env KDIR=$L/klem_expand_c_fam$j K=10 NP=10 MAXITS=6 RTAG=_c_fam$j bash $W/run_state_refine_toyc.sbatch 2>&1 | grep "k-means\|REFINED\|poses\|FAILED" | sed "s/^/[fam $j refine] /"
  RD=$(ls -td $L/state_refine_k10_c_fam$j/*_refine3D_states 2>/dev/null | head -1); [ -z "$RD" ] && { log "family $j: no refinement"; continue; }
  simple_exec prg=print_project_field projfile=$RD/p.simple oritype=ptcl3D mkdir=no 2>/dev/null > $L/refined_c_fam$j.txt
done
$PY - <<PYEOF
import numpy as np
rows = [l for l in open('$L/ptcl3D_r0.txt') if l.strip() and not l.startswith('#')]; N = len(rows); n = 0
for j in (1, 2):
    try: t = [l for l in open('$L/refined_c_fam%d.txt' % j) if l.strip() and not l.startswith('#')]
    except Exception: continue
    assert len(t) == N
    for i, l in enumerate(t):
        kv = dict(x.split('=', 1) for x in l.split() if '=' in x)
        if float(kv.get('state', 0)) > 0: rows[i] = l.replace('state=%s' % kv['state'], 'state=1.00000000'); n += 1
open('$RT', 'w').write(''.join(rows)); print('realigned table: %d rows, %d from the family refinements' % (N, n))
PYEOF
log "poses: $($PY $W/pose_error.py $RT $S/all.txt --tag realigned)"
log "4. children again at the refined poses (feedback recipe, balanced)"; env PARENT=$L/klem_junk ETAG=_cfb NFAM=2 NFREE=0 TABLE=$RT PMARG=${PM:-3} DEVW=0.05 PNOISE=1 BAL=${BAL:-0.5} LPS="20 20 20 14 14 14 10 10 10" SIGMAS="1.5 1.5 1.2 1.0 1.0 0.8 0.7 0.7 0.7" bash $W/run_expand_toyc.sbatch 2>&1 | grep "final:\|MERGED\|best-state\|Error\|Traceback" | sed 's/^/[feedback] /'
log "5. level 3 (balanced, focus, 12->8 A)"; env CB=$L/klem_expand_cfb_fam TABLE=$RT BAL=${BAL:-0.5} PMARG=${PM:-3} L3=_cl3 bash $W/run_tree_level3b_toyc.sbatch 2>&1 | grep "grandchild\|states: cross-half\|MERGED\|best-state\|Error\|Traceback" | sed 's/^/[level3] /'
log DONE

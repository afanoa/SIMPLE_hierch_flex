#!/bin/bash
#SBATCH -J rnasep_closed_subcls_members
#SBATCH -p csbdevel
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 03:00:00
#SBATCH --mail-type=NONE
#SBATCH -o /mnt/beegfs/elmlund/afan/rnasep_klem/%x.%j.out
# The user sees the closed state in the sub-class averages the 3D focal score calls most docked. Test: take every
# sub-class with shrunk mean focal score >= SHR and >= MINN members on the family-1 board ("closed sub-classes"), split
# their members into C = also in the docked top 30k by own score and B = NOT in it (the per-particle score missed them),
# reconstruct B, C and A = B+C at the family poses and score the contact. If B carries the contact, the per-particle
# score under-rates closed particles that the 2D averages catch. Montage of the chosen sub-classes with ids.
# env: FD (board dir), SHR (0.40), MINN (20), N1 (30000)
set -uo pipefail
export PATH=/mnt/nasapps/production/gcc/15.2.0/bin:$PATH LD_LIBRARY_PATH=/mnt/nasapps/production/gcc/15.2.0/lib64
for v in $(env | grep -o "^SIMPLE_[A-Z0-9_]*" | grep -v "^SIMPLE_PATH$\|^SIMPLE_QSYS"); do unset $v; done
B=$HOME/SIMPLE_kernel/snap/5e1e6aeee; export SIMPLE_PATH=$B PATH=$B/bin:$B/scripts:$PATH SIMPLE_QSYS=slurm SIMPLE_QSYS_PARTITION=csbdevel OMP_NUM_THREADS=8 TOY_SMPD=1.464
PY=/mnt/nasapps/production/python/3.12/bin/python3; A=/mnt/beegfs/elmlund/afan; W=$A/routes/toy55; P=$A/rnasep_klem/expand1/fam1_pipeline; KF=$P/klem_fb; FD=${FD:-$KF/closed_focal}; META=$A/rnasep_coupled/basis_span_all/meta_all.npz
CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs.mrc; [ -s $CAVGS ] || CAVGS=$A/rnasep_coupled/full692k_2d/1_cls_expansion/cls_expansion_cavgs_even.mrc
SHR=${SHR:-0.40}; MINN=${MINN:-20}; N1=${N1:-30000}; D=$FD/box_subcls/members; mkdir -p $D && cd $D; log(){ echo "[$(date '+%F %T') cm] $*"; }
$PY - <<PYEOF
import numpy as np, sys
sys.path.insert(0, '$W'); from model_class import read_mrc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
M = np.load('$META'); flags = M['flags'].astype(int); lab = M['lab692'].astype(int); N = len(flags); ref = np.where(flags == 1)[0]; pos = -np.ones(N, int); pos[ref] = np.arange(len(ref))
F = np.loadtxt('$FD/board/focal_scores.txt', comments='#'); fsc = F[:, 1]; R = np.load('$FD/board/klem_resp.npz'); r = R['r'].astype(float); pinds = R['pinds'].astype(int); ngrid = len(R['i1']); s = r[:, :ngrid] @ fsc[:ngrid]
sl = lab[pos[pinds - 1]]; top = np.zeros(N + 1, bool); top[pinds[np.argsort(-s)[:$N1]]] = True
T = np.loadtxt('$FD/box_subcls/subclass_focal_means.txt', comments='#'); u, cnt, mean, shrunk = T[:, 0].astype(int), T[:, 1].astype(int), T[:, 2], T[:, 3]
pick = u[(shrunk >= $SHR) & (cnt >= $MINN)]; print('closed sub-classes (shrunk mean >= %.2f, >= %d board members): %d of %d; members on the board %d' % ($SHR, $MINN, len(pick), len(u), int(np.isin(sl, pick).sum())))
inA = np.isin(sl, pick); inC = inA & top[pinds]; inB = inA & ~top[pinds]
print('  A = all their members %d; C = also in the docked top %d by own score %d (%.0f %%); B = missed by the per-particle score %d; B own-score median %+.2f vs C %+.2f' % (inA.sum(), $N1, inC.sum(), 100. * inC.sum() / max(inA.sum(), 1), inB.sum(), np.median(s[inB]) if inB.any() else float('nan'), np.median(s[inC]) if inC.any() else float('nan')))
st = np.zeros(N, int); st[pinds[inB] - 1] = 1; st[pinds[inC] - 1] = 2; np.savetxt('$D/states_BC.txt', st, fmt='%d'); sa = np.zeros(N, int); sa[pinds[inA] - 1] = 1; np.savetxt('$D/states_A.txt', sa, fmt='%d')
np.savetxt('$D/closed_subclasses.txt', np.c_[pick, cnt[np.isin(u, pick)], shrunk[np.isin(u, pick)]], fmt='%d %d %.3f', header='subclass n_on_board shrunk_mean_focal')
try:
    C = read_mrc('$CAVGS'); o = np.argsort(-shrunk[np.isin(u, pick)]); ids = pick[o][:48]; nr = int(np.ceil(len(ids) / 8.)); fig, ax = plt.subplots(nr, 8, figsize=(16, 2.1 * nr))
    for i, a in enumerate(np.atleast_1d(ax).ravel()):
        a.axis('off')
        if i < len(ids): k = ids[i]; a.imshow(C[k - 1], cmap='gray'); a.set_title('sub %d n%d %+.2f' % (k, cnt[u == k][0], shrunk[u == k][0]), fontsize=7)
    plt.suptitle('closed sub-classes by the 3D focal score (shrunk mean >= %.2f): average k = sub-class k' % $SHR); plt.tight_layout(); plt.savefig('$D/closed_subclasses_montage.png', dpi=110); print('montage -> $D/closed_subclasses_montage.png')
except Exception as e: print('montage skipped: %s' % e)
PYEOF
PS=$(ls -td $P/*_refine3D_states | head -1)/p.simple
rec(){ tag=$1; sf=$2; ns=$3; mkdir -p $D/rec_$tag && cd $D/rec_$tag && cp -f $PS p.simple; simple_exec prg=selection projfile=p.simple oritype=ptcl3D infile=$sf mkdir=no > sel.log 2>&1 || { log "selection FAILED ($tag)"; return 1; }
  simple_exec prg=reconstruct3D projfile=p.simple pgrp=c1 mskdiam=232 nstates=$ns nparts=6 nthr=8 > rec.log 2>&1 || { log "reconstruct3D FAILED ($tag): $(grep -a -m2 ERROR rec.log | cut -c1-200)"; return 1; }
  RD=$(ls -td ./*_reconstruct3D | head -1); mkdir -p maps; for w in 01 02; do for h in "" _even _odd; do cp -f $RD/recvol_state$w$h.mrc maps/ 2>/dev/null; done; done; cd $D; }
rec BC $D/states_BC.txt 2; rec A $D/states_A.txt 1
$PY $FD/lever/vis.py $FD/lever/maps/recvol_state01.mrc "docked top $N1 by own score" $D/rec_A/maps/recvol_state01.mrc "A: all members of the closed sub-classes" $D/rec_BC/maps/recvol_state02.mrc "C: members also in the top $N1" $D/rec_BC/maps/recvol_state01.mrc "B: members the per-particle score missed"
log DONE

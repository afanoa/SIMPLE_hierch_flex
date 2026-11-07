# pcg_focus invariant on a refine3D directory: every state map (merged, even, odd) equals the consensus background where the
# focus mask is zero, and differs from it only inside. usage: python3 focus_check.py <refine3D dir> <mask.mrc> <nstates>
import sys
import numpy as np, mrcfile

d, mfile, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
M = mrcfile.read(mfile).astype(np.float32)
out, ins = M <= 0., M > 0.5
for kind in ('', '_even', '_odd'):
    V = np.array([mrcfile.read('%s/recvol_state%02d%s.mrc' % (d, s, kind)).astype(np.float32) for s in range(1, k + 1)])
    scale = np.abs(V).max()
    dout = max(np.abs(V[i][out] - V[0][out]).max() for i in range(1, k)) / scale
    din = [np.sqrt(np.mean((V[i][ins] - V.mean(0)[ins]) ** 2)) / np.sqrt(np.mean(V.mean(0)[ins] ** 2)) for i in range(k)]
    print('maps%-5s: max |state - state 1| outside the mask / max|map| = %.3g; inside, rms(state - mean) / rms(mean) per state: %s'
          % (kind or '', dout, ' '.join('%.3f' % x for x in din)))
bg = mrcfile.read('%s/recvol_focus_bg.mrc' % d).astype(np.float32)
V = np.array([mrcfile.read('%s/recvol_state%02d.mrc' % (d, s)).astype(np.float32) for s in range(1, k + 1)])
print('merged maps vs recvol_focus_bg.mrc outside the mask: max |diff| / max|map| = %.3g' % (np.abs(V[:, out] - bg[out]).max() / np.abs(V).max()))

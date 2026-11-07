# pcg_focus invariant on a refine3D directory whose maps are at a cropped box: where the state maps agree exactly they
# must equal recvol_focus_bg.mrc, and the voxels where they differ must lie inside the focus mask (resampled to the crop)
# usage: python3 focus_check2.py <refine3D dir> <mask.mrc at the native box> <nstates>
import sys
import numpy as np, mrcfile
from scipy import ndimage

d, mfile, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
bg = mrcfile.read(d + '/recvol_focus_bg.mrc').astype(np.float32)
M = mrcfile.read(mfile).astype(np.float32)
Mz = ndimage.zoom(M, bg.shape[0] / M.shape[0], order=1)
for kind in ('', '_even', '_odd'):
    V = np.array([mrcfile.read('%s/recvol_state%02d%s.mrc' % (d, s, kind)).astype(np.float32) for s in range(1, k + 1)])
    same = np.all(V == V[0], axis=0)
    diff = ~same
    line = 'maps%-5s: states differ in %.2f %% of the box; of those voxels %.3f have resampled mask > 0.01, %.3f > 0.5' % (
        kind, 100. * diff.mean(), (Mz[diff] > 0.01).mean(), (Mz[diff] > 0.5).mean())
    line += '; mask > 0.5 voxels where states differ: %.3f' % diff[Mz > 0.5].mean()
    if kind == '':
        line += '; max |state - bg| where states agree / max|map| = %.3g' % (np.abs(V[:, same] - bg[same]).max() / np.abs(V).max())
    print(line)

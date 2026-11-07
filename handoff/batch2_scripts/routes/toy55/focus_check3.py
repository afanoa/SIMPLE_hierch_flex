# size of the state-to-background differences inside / outside the focus mask (resampled to the map box when cropped)
# usage: python3 focus_check3.py <refine3D dir> <mask.mrc at the native box> <nstates>
import sys
import numpy as np, mrcfile
from scipy import ndimage

d, mfile, k = sys.argv[1], sys.argv[2], int(sys.argv[3])
bg = mrcfile.read(d + '/recvol_focus_bg.mrc').astype(np.float32)
M = mrcfile.read(mfile).astype(np.float32)
Mz = M if M.shape == bg.shape else ndimage.zoom(M, bg.shape[0] / M.shape[0], order=1)
V = np.array([mrcfile.read('%s/recvol_state%02d.mrc' % (d, s)).astype(np.float32) for s in range(1, k + 1)])
D = np.abs(V - bg)
inside, outside, far = Mz > 0.5, Mz <= 0., ndimage.distance_transform_edt(Mz <= 0.01) > 5
print('box %d: |state - bg| max inside mask %.4g, outside (mask <= 0) %.4g, > 5 voxels from the mask %.4g; map max %.4g'
      % (bg.shape[0], D[:, inside].max(), D[:, outside].max() if outside.any() else 0., D[:, far].max() if far.any() else 0.,
         np.abs(V).max()))
print('   outside the mask: voxels where some state differs from bg %.2f %%, rms diff / rms inside %.3g'
      % (100. * (D[:, outside].max(0) > 0).mean() if outside.any() else 0.,
         np.sqrt((D[:, outside] ** 2).mean()) / np.sqrt((D[:, inside] ** 2).mean()) if outside.any() else 0.))

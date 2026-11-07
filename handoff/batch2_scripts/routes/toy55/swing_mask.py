# swept-domain focus mask of the swing30hs family: per-voxel std over the 31 swing conformer maps (conf200-230),
# smoothed at 4 A, thresholded at 0.2 of its maximum, largest component, dilated 6 A, 6 A cosine edge
# usage: python3 swing_mask.py <dir with conf2NN_swing_m160.mrc> <out.mrc> <region_5.mrc for the overlap report>
import glob, sys
import numpy as np, mrcfile
from scipy import ndimage

src, out, r5 = sys.argv[1], sys.argv[2], sys.argv[3]
fs = sorted(glob.glob(src + '/conf2[0-3][0-9]_swing_m160.mrc'))
V = np.array([mrcfile.read(f).astype(np.float32) for f in fs])
with mrcfile.open(fs[0]) as m:
    apix = float(m.voxel_size.x)
S = ndimage.gaussian_filter(V.std(0), 4.0 / apix)
B = S > 0.2 * S.max()
lab, n = ndimage.label(B)
if n > 1:
    B = lab == 1 + int(np.argmax(ndimage.sum(B, lab, range(1, n + 1))))
d = ndimage.distance_transform_edt(~B) * apix
M = np.where(d <= 6., 1., np.where(d < 12., 0.5 * (1. + np.cos(np.pi * (d - 6.) / 6.)), 0.)).astype(np.float32)
with mrcfile.new(out) as m:
    m.set_data(M)
    m.voxel_size = apix
R = mrcfile.read(r5).astype(np.float32) > .5
Mb = M > .5
print('%d swing maps (%s .. %s), apix %.3f, box %d' % (len(fs), fs[0].split('/')[-1], fs[-1].split('/')[-1], apix, M.shape[0]))
print('core %d voxels; mask > 0.5: %d voxels (%.2f %% of the box); region_5: %d voxels, %.2f inside the mask; mask inside region_5: %.2f'
      % (B.sum(), Mb.sum(), 100. * Mb.mean(), R.sum(), (R & Mb).sum() / max(1, R.sum()), (R & Mb).sum() / max(1, Mb.sum())))

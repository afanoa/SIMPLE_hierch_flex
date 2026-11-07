#!/usr/bin/env python3
"""land_cache.py: crop-C Fourier cache of a particle stack (lsamp_score2.py's build_cache without the real-space disc: independent pixel
noise for the pixel split). usage: land_cache.py STACKDIR OUT.npy [--crop 48]"""
import sys, argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('stack'); ap.add_argument('out'); ap.add_argument('--crop', type=int, default=48); a = ap.parse_args(); c = a.crop
h = np.frombuffer(open(a.stack + '/all.mrcs', 'rb').read(1024), np.int32); nx, nz, nsym = int(h[0]), int(h[2]), int(h[23])
I = np.memmap(a.stack + '/all.mrcs', np.float32, 'r', offset=1024 + nsym, shape=(nz, nx, nx)); sl = np.r_[0:c // 2, nx - c // 2:nx]
out = np.lib.format.open_memmap(a.out, 'w+', np.complex64, (nz, c, c))
for b in range(0, nz, 4000): blk = np.array(I[b:b + 4000], np.float32); out[b:b + len(blk)] = np.fft.fft2(np.fft.ifftshift(blk, axes=(1, 2)))[:, sl][:, :, sl]
out.flush(); print('cache %s: %d particles' % (a.out, nz))

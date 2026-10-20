#!/usr/bin/env python3
"""Deposited RNase P conformer maps (EMDB, 0.732 A, box 448/480) -> our lattice: Fourier-crop by 2 to 1.464 A, crop a
160 box around the density centroid, mirror (z-flip) into our hand. Docking to the consensus frame is done afterwards
with SIMPLE dock_volpair. usage: process_maps.py CSV EMDB_DIR OUT_DIR"""
import sys, os, gzip, csv, numpy as np
csvf, src, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)


def read_map_gz(fn):
    with gzip.open(fn, 'rb') as f:
        h = np.frombuffer(f.read(1024), np.int32); nx, ny, nz, mode, nsym = int(h[0]), int(h[1]), int(h[2]), int(h[3]), int(h[23])
        hf = h.view(np.float32); cella = hf[10:13]; f.read(nsym)
        dt = {2: np.float32, 0: np.int8, 1: np.int16, 6: np.uint16}[mode]
        a = np.frombuffer(f.read(nx * ny * nz * np.dtype(dt).itemsize), dt).reshape(nz, ny, nx).astype(np.float64)
    return a, float(cella[0]) / nx, (int(h[16]), int(h[17]), int(h[18]))


def write_mrc(fn, vol, smpd):
    v = np.ascontiguousarray(vol, np.float32); nz, ny, nx = v.shape
    h = np.zeros(256, np.int32); hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz); hf[10:13] = (nx * smpd, ny * smpd, nz * smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean()); h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little'); hf[54] = v.std()
    with open(fn, 'wb') as fo: fo.write(h.tobytes()); fo.write(v.tobytes())


def fourier_crop2(a):
    F = np.fft.fftshift(np.fft.fftn(a)); n = a.shape[0]; m = n // 2; o = (n - m) // 2
    return np.real(np.fft.ifftn(np.fft.ifftshift(F[o:o + m, o:o + m, o:o + m]))) / 8.


rows = [r for r in csv.DictReader(open(csvf))]
for r in rows:
    eid = r['emdb_id'].replace('EMD-', ''); fn = os.path.join(src, 'emd_%s.map.gz' % eid); of = os.path.join(out, 'conf%02d_emd%s_m160.mrc' % (int(r['conformer']), eid))
    if os.path.exists(of): continue
    a, px, axes = read_map_gz(fn); assert a.shape[0] == a.shape[1] == a.shape[2] and axes == (1, 2, 3), (a.shape, axes)
    b = fourier_crop2(a); smpd = px * 2
    # centroid of the positive density above 3 sd, crop 160 around it
    thr = b.mean() + 3 * b.std(); w = np.where(b > thr, b, 0.); zz, yy, xx = np.indices(b.shape); c = [int(round((w * q).sum() / w.sum())) for q in (zz, yy, xx)]
    B = 160; out3 = np.zeros((B, B, B)); 
    for ax in range(3): pass
    z0, y0, x0 = [ci - B // 2 for ci in c]
    zs, ys, xs = [slice(max(0, s), min(b.shape[0], s + B)) for s in (z0, y0, x0)]
    sub = b[zs, ys, xs]; oz, oy, ox = zs.start - z0, ys.start - y0, xs.start - x0
    out3[oz:oz + sub.shape[0], oy:oy + sub.shape[1], ox:ox + sub.shape[2]] = sub
    out3 = out3[::-1, :, :]          # mirror: our reconstructions are the other hand
    write_mrc(of, out3, smpd); print('%s: %s px %.3f -> %s @ %.3f, centroid %s, inside-box fraction %.3f' % (r['emdb_id'], a.shape, px, out3.shape, smpd, c, w[zs, ys, xs].sum() / w.sum()))

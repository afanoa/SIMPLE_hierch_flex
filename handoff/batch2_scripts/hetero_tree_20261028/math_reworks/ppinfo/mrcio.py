"""Minimal MRC I/O shared by the PfCRT micelle fixture scripts (volumes as v[z, y, x], stacks as s[n, y, x])."""
import numpy as np


def read_hdr(fn):
    h = np.fromfile(fn, np.int32, 256)
    return int(h[0]), int(h[1]), int(h[2]), int(h[23]), float(np.frombuffer(h[10:11].tobytes(), np.float32)[0]) / max(int(h[7]), 1)


def read_mrc(fn):
    nx, ny, nz, ns, smpd = read_hdr(fn)
    v = np.fromfile(fn, np.float32, offset=1024 + ns, count=nx * ny * nz).reshape(nz, ny, nx)
    return v.astype(np.float64), smpd


def stack(fn, mode='r'):
    nx, ny, nz, ns, _ = read_hdr(fn)
    return np.memmap(fn, np.float32, mode, offset=1024 + ns, shape=(nz, ny, nx))


def write_mrc(fn, v, smpd):
    v = np.ascontiguousarray(v, np.float32)
    nz, ny, nx = v.shape
    h = np.zeros(256, np.int32)
    hf = h.view(np.float32)
    h[0:3] = (nx, ny, nz); h[3] = 2; h[7:10] = (nx, ny, nz)
    hf[10:13] = (nx * smpd, ny * smpd, nz * smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    hf[19:22] = (v.min(), v.max(), v.mean())
    h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little')
    with open(fn, 'wb') as f:
        f.write(h.tobytes()); f.write(v.tobytes())


def write_stack_header(f, nx, ny, n, smpd):
    h = np.zeros(256, np.int32)
    hf = h.view(np.float32)
    h[0:3] = (nx, ny, n); h[3] = 2; h[7:10] = (nx, ny, n)
    hf[10:13] = (nx * smpd, ny * smpd, n * smpd); hf[13:16] = 90.; h[16:19] = (1, 2, 3)
    h[52] = int.from_bytes(b'MAP ', 'little'); h[53] = int.from_bytes(bytes([68, 65, 0, 0]), 'little')
    f.write(h.tobytes())


def fcrop3(v, n):
    """Fourier crop (or pad) a cubic volume to box n, preserving the real-space values' scale."""
    N = v.shape[0]
    F = np.fft.fftshift(np.fft.fftn(v))
    c, h = N // 2, n // 2
    if n <= N:
        G = F[c - h:c + h, c - h:c + h, c - h:c + h]
    else:
        G = np.zeros((n, n, n), complex); G[h - c:h + c, h - c:h + c, h - c:h + c] = F
    return np.real(np.fft.ifftn(np.fft.ifftshift(G))) * (n / N) ** 3


def lowpass3(v, smpd, lp, width=None):
    """Cosine-edged low-pass to lp Angstrom."""
    N = v.shape[0]
    f = np.fft.fftfreq(N, smpd)
    s = np.sqrt(f[:, None, None] ** 2 + f[None, :, None] ** 2 + f[None, None, :] ** 2)
    k = 1. / lp
    w = width if width is not None else 0.1 * k
    flt = np.clip(0.5 * (1 + np.cos(np.pi * np.clip((s - k) / w + 0.5, 0, 1))), 0, 1)
    flt[s < k - w / 2] = 1.
    return np.real(np.fft.ifftn(np.fft.fftn(v) * flt))


def read_pdb(fn, chains=None):
    xyz, info = [], []
    for ln in open(fn):
        if not ln.startswith(('ATOM', 'HETATM')):
            continue
        ch = ln[21]
        if chains is not None and ch not in chains:
            continue
        xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        info.append((ch, int(ln[22:26]), ln[12:16].strip(), ln[17:20]))
    return np.array(xyz), info

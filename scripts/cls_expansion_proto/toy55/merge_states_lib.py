import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from model_class import read_mrc
def leaf_ok(level1_dir, s, smpd=1.464, fscmax=float(__import__('os').environ.get('FSCMAX', '12')), minrec=int(__import__('os').environ.get('MINREC', '1000'))):
    """leaf support = the state's own half-map FSC (<= fscmax A) and its reconstructed count (>= minrec); junk reproduces across halves but does not reconstruct"""
    import os, numpy as np
    fe, fo = '%s/maps/recvol_state%02d_even.mrc' % (level1_dir, s), '%s/maps/recvol_state%02d_odd.mrc' % (level1_dir, s)
    if not (os.path.isfile(fe) and os.path.isfile(fo)): return True
    nrec = 0
    if os.path.isfile(level1_dir + '/clusters.txt'):
        for l in open(level1_dir + '/clusters.txt'):
            if not l.startswith('#') and int(l.split()[0]) == s: nrec = int(l.split()[3])
    if nrec and nrec < minrec: return False
    e, o = read_mrc(fe).astype(np.float32), read_mrc(fo).astype(np.float32); n = e.shape[0]; x = np.arange(n) - n // 2
    rr = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2); msk = np.clip((0.45 * n - rr) / (0.05 * n), 0, 1)
    E, O = np.fft.fftn(e * msk), np.fft.fftn(o * msk); f = np.fft.fftfreq(n); kx, ky, kz = np.meshgrid(f, f, f, indexing='ij'); sh = np.minimum(np.round(np.sqrt(kx**2 + ky**2 + kz**2) * n).astype(int), n // 2)
    num = np.bincount(sh.ravel(), np.real(np.conj(E) * O).ravel()); d1 = np.bincount(sh.ravel(), np.abs(E).ravel() ** 2); d2 = np.bincount(sh.ravel(), np.abs(O).ravel() ** 2)
    fsc = num / np.maximum(np.sqrt(d1 * d2), 1e-30); below = np.where(fsc[1:n // 2] < 0.143)[0]; sres = (below[0] + 1) if len(below) else n // 2 - 1
    return n * smpd / max(sres, 1) <= fscmax


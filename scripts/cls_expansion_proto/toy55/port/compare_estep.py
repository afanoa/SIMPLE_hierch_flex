#!/usr/bin/env python3
"""Compare the Fortran klem E-step stream output with the Python dump (klem_estep_ref.npz): the centred grid of cell 1
(even), the particle band coefficients F, the CTF C, the projection P of cell 1, and the log-likelihoods ll.
usage: compare_estep.py DUMP_DIR FORTRAN_OUT"""
import sys, glob, numpy as np
d, fo = sys.argv[1], sys.argv[2]; Z = np.load(d + '/klem_estep_ref.npz')
with open(fo, 'rb') as f:
    n, K, nband, nsamp, n3 = np.fromfile(f, np.int32, 5); ll = np.fromfile(f, np.float32, n * K).reshape(n, K)
    Fs = np.fromfile(f, np.complex64, nsamp * nband).reshape(nsamp, nband); Cs = np.fromfile(f, np.float32, nsamp * nband).reshape(nsamp, nband)
    Ps = np.fromfile(f, np.complex64, nsamp * nband).reshape(nsamp, nband); V = np.fromfile(f, np.complex64, n3 ** 3).reshape(n3, n3, n3)
print('fortran: n %d K %d nband %d nsamp %d n3 %d | python: ll %s nband %d n3 %d' % (n, K, nband, nsamp, n3, Z['ll'].shape, len(Z['kx']), int(Z['n3'])))
def rel(a, b): return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-30))
def cc(a, b): a = a.ravel(); b = b.ravel(); return float(np.real(np.vdot(a, b)) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-30))
# cell 1 even: python centred grid from the MRC it wrote
def read_mrc(p):
    h = np.fromfile(p, dtype=np.int32, count=256); nx, ny, nz = h[0], h[1], h[2]; return np.fromfile(p, dtype=np.float32, offset=1024 + h[23]).reshape(nz, ny, nx).astype(float)
m = read_mrc(d + '/cells_half/cell_01_even.mrc'); Vf = np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(m))); cb = m.shape[0] // 2; hb = n3 // 2; Vp = Vf[cb - hb:cb - hb + n3, cb - hb:cb - hb + n3, cb - hb:cb - hb + n3]
print('V cell 1 even: corr %.5f, rel diff %.4f (conj corr %.5f)' % (cc(V, Vp), rel(V, Vp), cc(np.conj(V), Vp)))
rv = np.median(np.abs(V[np.abs(Vp) > np.percentile(np.abs(Vp), 90)]) / np.abs(Vp[np.abs(Vp) > np.percentile(np.abs(Vp), 90)])); print('V scale fortran/python: median ratio %.6g (1/N3 = %.3g, 1/sqrt(N3) = %.3g, 1/N3 of the padded box %.3g)' % (rv, 1. / n3 ** 3, 1. / n3 ** 1.5, 1. / m.shape[0] ** 3))
ns = min(nsamp, len(Z['F200']) if 'F200' in Z else 0)
if ns:
    Fp, Cp = Z['F200'][:ns], Z['C200'][:ns]
    print('F (%d particles): corr %.5f, rel diff %.4f | conj corr %.5f | sign-flipped corr %.5f' % (ns, cc(Fs[:ns], Fp), rel(Fs[:ns], Fp), cc(np.conj(Fs[:ns]), Fp), cc(-Fs[:ns], Fp)))
    print('C: corr %.5f, rel diff %.4f (abs: %.5f)' % (cc(Cs[:ns], Cp), rel(Cs[:ns], Cp), cc(np.abs(Cs[:ns]), np.abs(Cp))))
    big = np.abs(Fp) > np.percentile(np.abs(Fp), 90); rf = np.median(np.abs(Fs[:ns][big]) / np.abs(Fp[big])); box = int(Z['box']); print('F scale fortran/python: median ratio %.6g (1/box^2 = %.3g, 1/box = %.3g)' % (rf, 1. / box ** 2, 1. / box))
    # ll recomputed from the fortran F/P with python-consistent scales: does the mismatch vanish?
    print('P cell 1 (fortran) vs F/(C scale) is not available from python; P stats: |P| mean fortran %.3g' % np.abs(Ps[:ns]).mean())
llp = Z['ll'][:n]; print('ll: corr over particles x cells %.5f; per-particle centred corr median %.4f; rel diff %.4f; argmax agreement %.3f' % (cc(ll, llp), np.median([cc(ll[i] - ll[i].mean(), llp[i] - llp[i].mean()) for i in range(min(n, 2000))]), rel(ll, llp), float(np.mean(ll.argmax(1) == llp.argmax(1)))))
sm = lambda x: np.exp(x - x.max(1, keepdims=True)) / np.exp(x - x.max(1, keepdims=True)).sum(1, keepdims=True); print('responsibilities (softmax of ll): mean |r_f - r_p| %.4f, argmax agreement %.3f' % (float(np.abs(sm(ll) - sm(llp)).mean()), float(np.mean(sm(ll).argmax(1) == sm(llp).argmax(1)))))

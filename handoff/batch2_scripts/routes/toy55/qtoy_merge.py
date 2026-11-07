#!/usr/bin/env python3
"""usage: qtoy_merge.py SIMDIR CSV TRUTH_OUT   (env SNR0 0.0012, BMAX 200, QSIG 0.35, CHUNK 4000)
qtoy merge: the noise-free per-class stacks (simulate_particles adds no noise at snr >= 5) concatenated in a seeded random order (11), each
particle given a quality draw independent of its state (seed 23: B ~ U(0, BMAX) A^2, SNR scale q ~ lognormal(0, QSIG) clipped to 0.4-2.5)
and one absolute noise level (seed 31): image = sqrt(q) x exp(-B s^2 / 4) x clean + N(0, sigma^2), sigma^2 = the mean clean-image variance of
the non-junk particles over the 2x zero-padded box / SNR0 (simulate_particles adds its snr noise to the padded image before clipping: the sites
toy's RNase P absolute noise). CSV columns: conformer, images, kind (arm | closed | holo), s1,
closed_arm. Writes SIMDIR/all.mrcs, SIMDIR/all.txt (truth poses + CTF, state = conformer, 0 junk) and TRUTH_OUT, one row per particle:
row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq (combo = conformer, A_domain = holo, C_closed = the
closed state: toy_score.py's simsites columns)."""
import os, sys, csv, numpy as np
S, CSVF, TO = sys.argv[1:4]; SNR0 = float(os.environ.get('SNR0', '0.0012')); BMAX = float(os.environ.get('BMAX', '200')); QSIG = float(os.environ.get('QSIG', '0.35'))
CH = int(os.environ.get('CHUNK', '4000')); SMPD = 1.464
R = list(csv.DictReader(open(CSVF))); names = ['conf%02d' % int(r['conformer']) for r in R] + ['blob', 'broken']
info = [(int(r['conformer']), r['kind'], float(r['s1']), float(r['closed_arm']), 0) for r in R] + [(0, 'junk', np.nan, np.nan, 1), (0, 'junk', np.nan, np.nan, 2)]
def kv(line): return dict(t.split('=', 1) for t in line.split() if '=' in t)
mm, oris, cls = [], [], []
for k, nm in enumerate(names):
    fn = '%s/stk_%s.mrcs' % (S, nm); h = np.fromfile(fn, dtype=np.int32, count=256); nx, ny, nz, nsym = (int(x) for x in (h[0], h[1], h[2], h[23]))
    mm.append(np.memmap(fn, dtype=np.float32, mode='r', offset=1024 + nsym, shape=(nz, ny, nx)))
    o = [x.rstrip('\n') for x in open('%s/oris_%s.txt' % (S, nm)) if x.strip()]; assert len(o) == nz, (nm, nz, len(o)); oris += o; cls += [(k, i) for i in range(nz)]
cls = np.array(cls); N = len(cls); perm = np.random.default_rng(11).permutation(N); src = cls[perm]
rq = np.random.default_rng(23); B = rq.uniform(0., BMAX, N); Q = np.clip(np.exp(rq.normal(0., QSIG, N)), 0.4, 2.5); rn = np.random.default_rng(31)
gen = np.where([info[k][4] == 0 for k in src[:, 0]])[0]; smp = np.random.default_rng(5).choice(gen, min(4000, len(gen)), replace=False)
def vpad(x): x = x.astype(np.float64); P = 4. * x.size; return (x * x).sum() / P - (x.sum() / P) ** 2
cv = np.array([vpad(mm[src[j, 0]][src[j, 1]]) for j in smp]); sig = np.sqrt(cv.mean() / SNR0)
yy, xx = np.mgrid[:ny, :nx]; rr = np.hypot(xx - nx / 2, yy - ny / 2); cor = np.array([mm[src[j, 0]][src[j, 1]][rr > 0.45 * nx].var() / mm[src[j, 0]][src[j, 1]].var() for j in smp[:200]])
print('noise-free check: corner / whole-image variance of the clean stacks, median %.4f (an snr-10 noise would give about 0.1)' % np.median(cor))
print('absolute noise sigma %.4f = sqrt(mean clean variance %.5f / SNR0 %.4f); B ~ U(0, %g), q median %.2f (P2.5 %.2f, P97.5 %.2f)' % (sig, cv.mean(), SNR0, BMAX, np.median(Q), *np.percentile(Q, (2.5, 97.5))))
s2 = np.fft.fftfreq(ny, SMPD)[:, None] ** 2 + np.fft.rfftfreq(nx, SMPD)[None, :] ** 2
hdr = np.fromfile('%s/stk_%s.mrcs' % (S, names[0]), dtype=np.int32, count=256).copy(); hdr[2] = N; hdr[9] = N; hdr[23] = 0
with open(S + '/all.mrcs', 'wb') as f:
    f.write(hdr.tobytes())
    for c0 in range(0, N, CH):
        js = np.arange(c0, min(N, c0 + CH)); clean = np.empty((len(js), ny, nx), np.float32)
        for k in np.unique(src[js, 0]):
            w = np.where(src[js, 0] == k)[0]; loc = src[js[w], 1]; o = np.argsort(loc); clean[w[o]] = mm[k][loc[o]]
        F = np.fft.rfft2(clean) * (np.sqrt(Q[js])[:, None, None] * np.exp(-B[js][:, None, None] * s2[None] / 4.))
        f.write((np.fft.irfft2(F, s=(ny, nx)) + rn.normal(0., sig, clean.shape)).astype(np.float32).tobytes())
off = np.r_[0, np.cumsum([m.shape[0] for m in mm])]
with open(S + '/all.txt', 'w') as fa, open(TO, 'w') as ft:
    ft.write('# row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq\n')
    for j in range(N):
        k, i = src[j]; conf, kind, s1, ca, jt = info[k]; d = kv(oris[off[k] + i]); d['state'] = str(conf)
        fa.write(' '.join('%s=%s' % (a, b) for a, b in d.items()) + '\n')
        ft.write('%d %d %d %d 0 %d %d %d %s %s %d %.2f %.4f\n' % (j + 1, conf, int(jt > 0), int(kind == 'holo'), int(kind == 'closed'), jt, conf, '%.2f' % s1, '%.3f' % ca, int(kind == 'holo'), B[j], Q[j]))
cnt = {}; [cnt.__setitem__(names[k], cnt.get(names[k], 0) + 1) for k in src[:, 0]]
print('merged %d particles: %s' % (N, ', '.join('%s %d' % (k, v) for k, v in cnt.items())))

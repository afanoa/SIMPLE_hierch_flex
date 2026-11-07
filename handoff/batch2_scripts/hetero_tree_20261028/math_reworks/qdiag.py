import numpy as np, sys
COLS = 'row combo junk A_domain B_blob C_closed jtype conformer s1 closed_arm holo bfac snrq'.split()
for name, truth, d in [('qtoynull root', '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoynull/truth_qtoy.txt', 'flag2/base_qtoynull_r'),
                       ('qtoy2 root', '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt', 'flag2/base_qtoy2_r'),
                       ('qtoy2 r1', '/mnt/beegfs/elmlund/afan/routes/toy55/simqtoy2/truth_qtoy.txt', 'flag2/base_qtoy2_r1')]:
    t = np.loadtxt(truth, ndmin=2); C = {c: t[:, COLS.index(c)] for c in COLS}; rows = C['row'].astype(int)
    ch = np.loadtxt(d + '/split_children_pre.txt', dtype=int)[rows - 1]
    w = np.loadtxt(d + '/ptcl_weights.txt', comments='#'); wr = w[:, 0].astype(int); wq = {int(r): (a, b, c, e) for r, a, b, c, e in zip(wr, w[:, 1], w[:, 2], w[:, 3], w[:, 4])}
    corr = np.array([wq.get(r, (0, 0, 0, np.nan))[3] for r in rows]); lnv = np.array([wq.get(r, (0, 0, np.nan, 0))[2] for r in rows]); con = np.array([wq.get(r, (0, np.nan, 0, 0))[1] for r in rows])
    gen = (C['junk'] == 0) & (ch > 0)
    print('==', name, '| groups', [int((ch == g).sum()) for g in range(1, ch.max() + 1)])
    for q in ['snrq', 'bfac', 'junk']:
        x = C[q]; sw = np.sqrt(np.mean([x[(ch == g) & (gen if q != 'junk' else ch > 0)].var() for g in range(1, ch.max() + 1)]))
        m = x[ch > 0].mean() if q == 'junk' else x[gen].mean()
        print('   %-5s group means %s | d %s (within sd %.4f)' % (q, ' '.join('%.4f' % x[(ch == g) & (gen if q != 'junk' else ch > 0)].mean() for g in range(1, ch.max() + 1)),
              ' '.join('%+.2f' % ((x[(ch == g) & (gen if q != 'junk' else ch > 0)].mean() - m) / sw) for g in range(1, ch.max() + 1)), sw))
    ok = gen & np.isfinite(corr)
    print('   corr(snrq, stage-B corr) %.3f  corr(snrq, log noise var) %.3f  corr(snrq, contrast) %.3f  corr(bfac, stage-B corr) %.3f' % (
        np.corrcoef(C['snrq'][ok], corr[ok])[0, 1], np.corrcoef(C['snrq'][ok], lnv[ok])[0, 1], np.corrcoef(C['snrq'][ok], con[ok])[0, 1], np.corrcoef(C['bfac'][ok], corr[ok])[0, 1]))
    # view / pose association: the split test strata are view x defocus; use s1 (a truth column) as a crude proxy of pose if present
    print('   s1 group means %s' % ' '.join('%.3f' % C['s1'][(ch == g) & gen].mean() for g in range(1, ch.max() + 1)))

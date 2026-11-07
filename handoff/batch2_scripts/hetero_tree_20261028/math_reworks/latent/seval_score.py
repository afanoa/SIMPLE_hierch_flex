#!/usr/bin/env python3
"""seval_score.py: a refine=state_eval run against the stage-B labels it started from and against the truth: the label agreement, and
for both labelings the truth readout of latent_states_v4 (qtoy: NMI with the combo and each truth class's best-state F1; closed: the
closed reference rows' best enrichment).
usage: seval_score.py NODE BEFORE_DUMP AFTER_DUMP TRUTHSPEC"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS


def states(fn):
    return np.array([int(float(dict(t.split('=', 1) for t in ln.split() if '=' in t).get('state', '0')))
                     for ln in open(fn) if 'e1=' in ln])


node, fb, fa, spec = sys.argv[1:5]
sb, sa = states(fb), states(fa)
didx, prow = LS.node_select(len(sb), node)
sb, sa = sb[didx], sa[didx]
act = sb > 0
print('label agreement with stage B %.4f over %d active particles' % ((sa[act] == sb[act]).mean(), act.sum()))
kind, T = LS.load_truth(spec, node, prow)
T = {k: v[act] for k, v in T.items()}
for tag, s in (('stage B   ', sb[act]), ('state_eval', sa[act])):
    rd = LS.readout(kind, T, s)
    line = ' '.join('%s %.4f' % (k, rd[k]) for k in ('nmi_combo', 'ari_combo', 'best_enrich') if k in rd)
    if kind == 'qtoy':
        f1 = []
        for c in np.unique(T['combo']):
            mc = T['combo'] == c
            if mc.mean() < 0.01:
                continue
            best = max(2 * (mc & (s == u)).sum() / (mc.sum() + (s == u).sum()) for u in np.unique(s))
            f1.append('%d:%.3f' % (c, best))
        line += '  best-state F1 ' + ' '.join(f1)
    print('%s  %s' % (tag, line))

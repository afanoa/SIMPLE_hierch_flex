#!/usr/bin/env python3
"""f2r_ab_readout.py: the stage-B A/B's labels against the truth and each other. Each arm is a run directory with b_ptcl3D.txt
(the stage-B output's ptcl3D dump) and split_children.txt (the split step's verdict); the production stage B (NODE/stageB) can be
given as an arm too (its dump under --dumps/<key>/b_ptcl3D.txt from the f2r runs).
usage: f2r_ab_readout.py --node NODE --truth KIND:ARGS --arm NAME=DIR_OR_DUMP [--arm ...]"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import latent_states_v4 as LS
from ppinfo_site import read_dump


def best_f1(groups, truth, v):
    """the best group for truth value v: F1, recall and purity (the positive-control toy's closed / holo readout)."""
    pos = truth == v
    best = (0., 0., 0.)
    for u in np.unique(groups[groups > 0]):
        m = groups == u
        tp = float((m & pos).sum())
        if tp == 0:
            continue
        rec, pur = tp / pos.sum(), tp / m.sum()
        f1 = 2 * rec * pur / (rec + pur)
        if f1 > best[0]:
            best = (f1, rec, pur)
    return best


def kind_line(groups, T):
    out = ''
    for v, name in ((18, 'closed'), (101, 'holo')):
        if (T['combo'] == v).any():
            out += ' %s F1 %.2f (recall %.2f, purity %.2f)' % ((name,) + best_f1(groups, T['combo'], v))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--node'); ap.add_argument('--truth', default='none'); ap.add_argument('--arm', action='append', default=[])
    a = ap.parse_args()
    rows = np.loadtxt(os.path.join(a.node, 'rows.txt'), dtype=int, ndmin=1)
    kind, T = LS.load_truth(a.truth, a.node, rows)
    labs = {}
    for arm in a.arm:
        name, path = arm.split('=', 1)
        fn = path if path.endswith('.txt') else os.path.join(path, 'b_ptcl3D.txt')
        d = read_dump(fn)
        st = d['state'].astype(int)
        if len(st) > len(rows) and rows.max() <= len(st):
            st = st[rows - 1]
        if len(st) != len(rows):
            raise SystemExit('%s: %d rows for a node of %d' % (fn, len(st), len(rows)))
        labs[name] = st
        pops = ' '.join('%d:%d' % (u, (st == u).sum()) for u in np.unique(st))
        line = '%-10s states %s' % (name, pops)
        if kind == 'closed':
            rd = LS.readout(kind, T, st)
            line += ' | closed enrichment of the best state %.3f (groups %s, separates %s)' % (rd['best_enrich'], rd['groups'], rd['separates'])
        elif kind == 'qtoy':
            rd = LS.readout(kind, T, st)
            line += ' | NMI %.3f ARI %.3f' % (rd['nmi_combo'], rd['ari_combo']) + kind_line(st, T)
        sc = os.path.join(path, 'split_children.txt') if not path.endswith('.txt') else ''
        if sc and os.path.exists(sc):
            ch = np.loadtxt(sc, dtype=int, ndmin=1)
            if len(ch) > len(rows) and rows.max() <= len(ch):
                ch = ch[rows - 1]
            line += ' | split children %s' % ' '.join('%d:%d' % (u, (ch == u).sum()) for u in np.unique(ch))
            if kind == 'qtoy' and len(ch) == len(rows) and len(np.unique(ch[ch > 0])) > 1:
                line += ' (children NMI %.3f%s)' % (LS.nmi(ch[ch > 0], T['combo'][ch > 0]), kind_line(ch, T))
        print(line)
    names = list(labs)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            x, y = labs[names[i]], labs[names[j]]
            m = (x > 0) & (y > 0)
            print('%s vs %s: same label %.3f of %d assigned (NMI %.3f)' % (names[i], names[j], (x[m] == y[m]).mean(), int(m.sum()), LS.nmi(x[m], y[m])))


if __name__ == '__main__':
    main()

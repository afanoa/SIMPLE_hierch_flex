#!/usr/bin/env python3
"""site_split gate 1: the Fortran worker's per-site calls and labels against the Python gate replays.

usage: xs_gate1_compare.py LIST RUNDIR GATESX [TABLE]
LIST    the xs_gate1.list lines "KEY NODE MSK PG"; RUNDIR holds KEY/summary.txt and KEY/sN/labels.txt (all project rows, 0 = not in
        the node); GATESX/runs/<KEY>_s0N (or _i0N)/labels_final.txt holds the Python labels in the node's rows.txt order.
TABLE   optional "KEY site gapw call" lines from the Python table, for the GAPw-within-.05 check.
"""
import os, sys
import numpy as np

def fnum(x):
    try: return float(x)
    except ValueError: return float('nan')

def read_summary(fn):
    out = {}
    if not os.path.exists(fn): return out
    for line in open(fn):
        f = line.split()
        if line.startswith('#') or len(f) < 13: continue
        out[f[0]] = dict(n1=int(f[1]), n2=int(f[2]), z1=fnum(f[3]), r1=fnum(f[4]), lsep=fnum(f[5]), ang=fnum(f[6]),
                         nb=int(f[7]), gapw=fnum(f[8]), gape=fnum(f[9]), gapo=fnum(f[10]), call=f[11], fails=f[12])
    return out

def agreement(a, b):
    """fraction agreeing under the better of the two label pairings (2 groups)"""
    same = np.mean(a == b); swap = np.mean(a == (3 - b))
    return max(same, swap)

def main():
    lst, rundir, gatesx = sys.argv[1:4]
    table = {}
    if len(sys.argv) > 4:
        for line in open(sys.argv[4]):
            f = line.split()
            if line.startswith('#') or len(f) < 4: continue
            table[(f[0], f[1])] = (float(f[2]), f[3], f[4] if len(f) > 4 else '-', f[5] if len(f) > 5 else '-',
                                   f[6] if len(f) > 6 else '-', f[8] if len(f) > 8 else '-')
    print(f"{'key':12s} {'site':5s} {'n1':>6s} {'n2':>6s} {'z1':>7s} {'r1':>5s} {'lsep':>6s} {'ang':>5s} {'nb':>5s} {'gapw':>6s} {'call':13s} "
          f"{'fails':12s} {'py_gapw':>7s} {'py_call':9s} {'py_z1':>5s} {'py_r1':>5s} {'py_L':>4s} {'py_B':>4s} {'agree':>6s}")
    ndisc = 0; nsites = 0
    for line in open(lst):
        f = line.split()
        if len(f) < 4: continue
        key, node = f[0], f[1]
        summ = read_summary(os.path.join(rundir, key, 'summary.txt'))
        rows = None
        rf = os.path.join(node, 'rows.txt')
        if os.path.exists(rf):
            rows = np.loadtxt(rf, dtype=int, ndmin=1)
        for tag, s in summ.items():
            nsites += 1
            if s['call'] == 'DISCRETE': ndisc += 1
            k = int(tag[1:])
            pyk = f"{key}_{'s' if tag[0] == 's' else 'i'}{k:02d}"
            pyd = os.path.join(gatesx, 'runs', pyk)
            if not os.path.isdir(pyd):
                alt = os.path.join(gatesx, 'runs', f"{key}_i{k:02d}")
                pyd = alt if os.path.isdir(alt) else pyd
            if key == 'igmfix_root' and k == 6: pyd = os.path.join(gatesx, 'runs', 'v_fix_a1')
            agree = float('nan')
            lf = os.path.join(rundir, key, tag, 'labels.txt'); pf = os.path.join(pyd, 'labels_final.txt')
            if os.path.exists(lf) and os.path.exists(pf) and rows is not None:
                lab = np.loadtxt(lf, dtype=int, ndmin=1); py = np.loadtxt(pf, dtype=int, ndmin=1)
                if lab.max() > 0:
                    nz = np.flatnonzero(lab) + 1
                    if len(nz) == len(py) and len(rows) == len(py):
                        a = lab[rows - 1]
                        agree = agreement(a, py) if np.all(a > 0) else agreement(lab[nz - 1], py)
                    elif len(nz) == len(py):
                        agree = agreement(lab[nz - 1], py)
            pg, pc, pz, pr, pl, pb = table.get((key, tag), (float('nan'), '-', '-', '-', '-', '-'))
            print(f"{key:12s} {tag:5s} {s['n1']:6d} {s['n2']:6d} {s['z1']:7.2f} {s['r1']:5.2f} {s['lsep']:6.1f} {s['ang']:5.1f} {s['nb']:5d} "
                  f"{s['gapw']:6.3f} {s['call']:13s} {s['fails']:12s} {pg:7.3f} {pc:9s} {pz:>5s} {pr:>5s} {pl:>4s} {pb:>4s} {agree:6.3f}")
    print(f"\n{nsites} sites, {ndisc} DISCRETE")

if __name__ == '__main__':
    main()

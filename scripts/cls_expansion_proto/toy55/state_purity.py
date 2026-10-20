#!/usr/bin/env python3
"""Per-conformer purity / recall of the state labels in a SIMPLE ptcl3D table (print_project_field) against the toy
truth; Hungarian one-state-per-conformer purity. usage: state_purity.py ptcl3D.txt gt_particles.txt [--tag T]"""
import sys, argparse, numpy as np
from scipy.optimize import linear_sum_assignment
ap = argparse.ArgumentParser(); ap.add_argument('table'); ap.add_argument('gt'); ap.add_argument('--tag', default=''); a = ap.parse_args()
gt = {}
for l in open(a.gt):
    if l.strip() and not l.startswith('#'): w = l.split(); gt[int(w[0])] = int(w[1])
st = []
for l in open(a.table):
    kv = dict(t.split('=', 1) for t in l.split() if '=' in t); st.append(int(float(kv.get('state', 0))))
st = np.array(st); conf = np.array([gt.get(i + 1, 0) for i in range(len(st))]); m = (st > 0) & (conf > 0); cs = sorted(set(conf[m])); ss = sorted(set(st[m]))
C = np.zeros((len(ss), len(cs)))
for i, s_ in enumerate(ss):
    for j, c in enumerate(cs): C[i, j] = np.sum((st == s_) & (conf == c) & m)
P = C / np.maximum(C.sum(1, keepdims=True), 1); R = C / np.maximum(C.sum(0, keepdims=True), 1); base = C.sum(0) / C.sum()
r_, c_ = linear_sum_assignment(-P); pur = np.zeros(len(cs)); rec = np.zeros(len(cs)); pur[c_] = P[r_, c_]; rec[c_] = R[r_, c_]
print('%s: %d states x %d conformers, %d particles; per conformer (one state each, Hungarian) purity/recall: %s | mean purity %.2f (min %.2f), mean recall %.2f; conformers with a state whose majority they are: %d/%d' % (
    a.tag, len(ss), len(cs), int(C.sum()), ' '.join('%d:%.2f/%.2f' % (c, pur[j], rec[j]) for j, c in enumerate(cs)), pur.mean(), pur.min(), rec.mean(), len(set(P.argmax(1))), len(cs)))

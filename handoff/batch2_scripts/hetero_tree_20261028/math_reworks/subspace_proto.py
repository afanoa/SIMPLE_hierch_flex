#!/usr/bin/env python3
"""Prototype of the audit's 6.3 endpoint test on site_rank's dumped hold-out-band Gram matrices.

For every site dir given: read site_gram.txt (even.even, odd.odd, even.odd Grams of the endpoint deviation maps and of
each final state's seed-set deviation), compute the single-fit cross-half tests (self reproducibility, residual
reproducibility given another direction = principal-angle test, disattenuated correlation), group the endpoints into
poles, and print them beside the halves verdicts recorded in site_halves.txt / site_pole_rule.txt / site_poles.txt.
"""
import sys, math, os, re
import numpy as np

KSIG = 3.0

def read_gram(path):
    with open(path) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    lines = [l for l in lines if not l.startswith('#')]
    neff, klo, khi = lines[0].split()[:3]
    neff = float(neff)
    names = lines[1].split()
    m = len(names)
    rows = [list(map(float, l.split())) for l in lines[2:2 + 3 * m]]
    gee = np.array(rows[:m]); goo = np.array(rows[m:2 * m]); geo = np.array(rows[2 * m:3 * m])
    return neff, (int(klo), int(khi)), names, gee, goo, geo

def fisher_z(r, n):
    r = max(min(r, 0.999999), -0.999999)
    return math.atanh(r) * math.sqrt(max(n - 3.0, 1.0))

class Grams:
    def __init__(self, gee, goo, geo, neff):
        self.gee, self.goo, self.geo, self.neff = gee, goo, geo, neff
        self.m = gee.shape[0]
    def ip(self, half1, u, half2, v):
        """inner product of combination u (weights over maps) in half1 with v in half2"""
        if half1 == 'e' and half2 == 'e': G = self.gee
        elif half1 == 'o' and half2 == 'o': G = self.goo
        elif half1 == 'e' and half2 == 'o': G = self.geo
        else: G = self.geo.T
        return float(u @ G @ v)
    def unit(self, a):
        """weight vector of map a normalised to unit norm in each half"""
        ue = np.zeros(self.m); ue[a] = 1.0 / math.sqrt(self.gee[a, a])
        uo = np.zeros(self.m); uo[a] = 1.0 / math.sqrt(self.goo[a, a])
        return ue, uo
    def group_dir(self, members):
        """mean of the unit-normalised member deviations, per half"""
        ue = np.zeros(self.m); uo = np.zeros(self.m)
        for s in members:
            e, o = self.unit(s); ue += e / len(members); uo += o / len(members)
        return ue, uo
    def cross_corr(self, a, b):
        return self.geo[a, b] / math.sqrt(self.gee[a, a] * self.goo[b, b])
    def self_z(self, a):
        return fisher_z(self.cross_corr(a, a), self.neff)
    def dirs(self, groups):
        """unit-mean direction of each group (list of member lists), per half: arrays (m, k)"""
        k = len(groups); We = np.zeros((self.m, k)); Wo = np.zeros((self.m, k))
        for j, g in enumerate(groups):
            ue, uo = self.group_dir(g); We[:, j] = ue; Wo[:, j] = uo
        return We, Wo
    def residual(self, b, groups):
        """cross-half correlation of b's residual after least-squares projection onto the span of the groups' directions,
        done in each half with that half's directions"""
        if isinstance(groups[0], (int, np.integer)): groups = [[g] for g in groups]
        We, Wo = self.dirs(groups)
        be = np.zeros(self.m); be[b] = 1.0; bo = be.copy()
        Ge = We.T @ self.gee @ We; Go = Wo.T @ self.goo @ Wo
        ce = np.linalg.solve(Ge, We.T @ self.gee @ be); co = np.linalg.solve(Go, Wo.T @ self.goo @ bo)
        re_ = be - We @ ce; ro_ = bo - Wo @ co
        num = self.ip('e', re_, 'o', ro_)
        den = math.sqrt(max(self.ip('e', re_, 'e', re_), 1e-300) * max(self.ip('o', ro_, 'o', ro_), 1e-300))
        r = num / den
        frac = self.ip('e', re_, 'e', re_) / self.gee[b, b]
        return r, fisher_z(r, self.neff - len(groups)), frac
    def disatt(self, a, b):
        return 0.5 * (self.cross_corr(a, b) + self.cross_corr(b, a)) / math.sqrt(max(self.cross_corr(a, a) * self.cross_corr(b, b), 1e-12))

def old_verdicts(d):
    """endpoint -> (pole, accept z row, accept z col); merges (a,b) yes; pole rule rows"""
    ends = {}; merges = {}; rule = []
    p = os.path.join(d, 'site_halves.txt')
    if os.path.exists(p):
        for l in open(p):
            w = l.split()
            if not w or l.startswith('#'): continue
            if w[0] == 'accept':
                ends.setdefault(w[1], {})['z' + w[2]] = float(w[-2])
            elif w[0] == 'merge':
                merges[(w[1], w[2])] = (w[-1] == 'yes', float(w[-2]))
            elif '|' in l and w[0] not in ('quarter', 'twin', 'pooled'):
                ends.setdefault(w[0], {})['pole'] = int(w[-1])
                ends[w[0]]['rel'] = (float(w[3]), float(w[4]))
    p = os.path.join(d, 'site_pole_rule.txt')
    if os.path.exists(p):
        for l in open(p):
            w = l.split()
            if not w or l.startswith('#'): continue
            rule.append(dict(pole=int(w[0]), seeds=int(w[1]), share=float(w[2]), zr=float(w[3]), zc=float(w[4]), fin=int(w[5]), kept=w[-1]))
    poles = [l.strip() for l in open(os.path.join(d, 'site_poles.txt'))] if os.path.exists(os.path.join(d, 'site_poles.txt')) else []
    starts = {}
    p = os.path.join(d, 'site_starts.txt')
    if os.path.exists(p):
        for l in open(p):
            w = l.split()
            if not w or l.startswith('#'): continue
            starts[w[0]] = dict(nseed=int(w[1]), nstar=int(w[2]), pole=int(w[3]), share=float(w[5]), hpk=float(w[6]), se=float(w[7]), epm=int(w[8]))
    return ends, merges, rule, poles, starts

def analyse(d, ksig=KSIG, verbose=True):
    neff, band, names, gee, goo, geo = read_gram(os.path.join(d, 'site_gram.txt'))
    G = Grams(gee, goo, geo, neff)
    ends, merges, rule, poles, starts = old_verdicts(d)
    iend = [i for i, n in enumerate(names) if not n.startswith('state')]
    ist  = [i for i, n in enumerate(names) if n.startswith('state')]
    out = []
    out.append(f'## {d}   n_eff {neff:.1f}  band {band[0]}..{band[1]}  (old) {" ; ".join(poles)}')
    # 1. endpoint self reproducibility
    acc = {}
    for a in iend:
        r = G.cross_corr(a, a); z = G.self_z(a); acc[a] = z >= ksig
        o = ends.get(names[a], {})
        out.append(f'  end {names[a]:5s} r_eo {r:6.3f} z_self {z:6.2f} {"ACC" if acc[a] else "rej"}'
                   f' | old pole {o.get("pole","-")} accept z row/col {o.get("zrow",float("nan")):6.2f}/{o.get("zcol",float("nan")):6.2f}'
                   f' hpk {starts.get(names[a],{}).get("hpk",float("nan")):.3f}')
    # 2. pairwise: disattenuated correlation and residual z both ways
    cand = [a for a in iend if acc[a]]
    pair = {}
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            rho = G.disatt(a, b)
            rb, zb, fb = G.residual(b, [a]); ra, za, fa = G.residual(a, [b])
            same = rho > 0 and zb < ksig and za < ksig   # anti-correlated ends are distinct poles by construction
            pair[(a, b)] = (rho, za, zb, same)
            om = merges.get((names[a], names[b]), merges.get((names[b], names[a]), (None, float('nan'))))
            out.append(f'  pair {names[a]:5s} {names[b]:5s} rho_hat {rho:6.3f} z_res(b|a) {zb:6.2f} z_res(a|b) {za:6.2f}'
                       f' resid share {fb:.2f}/{fa:.2f} -> {"SAME" if same else "distinct"} | old merge {om[0]} ({om[1]:.2f})')
    # 3. greedy grouping, closest first, complete linkage on the pair verdicts
    grp = {a: a for a in cand}
    order = sorted(pair.items(), key=lambda kv: -kv[1][0])
    for (a, b), (rho, za, zb, same) in order:
        if grp[a] == grp[b] or not same: continue
        A = [x for x in cand if grp[x] == grp[a]]; B = [x for x in cand if grp[x] == grp[b]]
        if all(pair.get((min(x, y), max(x, y)), (0, 0, 0, False))[3] for x in A for y in B):
            for x in B: grp[x] = grp[a]
    groups = {}
    for a in cand: groups.setdefault(grp[a], []).append(a)
    newpoles = list(groups.values())
    out.append('  new poles: ' + ' ; '.join('+'.join(names[a] for a in g) for g in newpoles) if newpoles else '  new poles: none')
    # 3b. redundancy: a pole whose deviation is explained by the span of the other positively correlated poles is dropped
    def pole_dir_corr(g, h):
        return np.mean([G.disatt(a, b) for a in g for b in h])
    keep = list(range(len(newpoles)))
    for i, g in sorted(enumerate(newpoles), key=lambda ig: np.mean([G.self_z(a) for a in ig[1]])):
        others = [newpoles[j] for j in keep if j != i and pole_dir_corr(g, newpoles[j]) > 0]
        if len(others) >= 2:
            zr = min(G.residual(a, others)[1] for a in g)
            out.append(f'  redundancy: pole {"+".join(names[a] for a in g)} vs span of {len(others)} positive poles: min z_res {zr:5.2f} -> {"DROP" if zr < ksig else "kept"}')
            if zr < ksig: keep.remove(i)
    newpoles = [newpoles[j] for j in keep]
    # 4. the states (seed-set deviations scored without EM): self z and distinctness vs the other positively correlated states
    for s in ist:
        r = G.cross_corr(s, s); z = G.self_z(s)
        others = [t for t in ist if t != s and G.disatt(s, t) > 0]
        zd = min([max(G.residual(s, [t])[1], G.residual(t, [s])[1]) for t in others]) if others else float('inf')
        m = re.match(r'state(\d+)', names[s]); k = int(m.group(1))
        o = next((x for x in rule if x['fin'] == k), None)
        neg = [names[t] for t in ist if t != s and G.disatt(s, t) <= 0]
        out.append(f'  {names[s]:28s} r_eo {r:6.3f} z_self {z:6.2f} z_dist {zd:6.2f} (anti/orth: {",".join(neg) or "-"}) | old rule: ' +
                   (f'share {o["share"]:.3f} z {o["zr"]:.2f}/{o["zc"]:.2f} {o["kept"]}' if o else 'rest/none'))
    # 5. rank of the accepted endpoint set: eigenvalues of the symmetrised cross-half correlation against a null
    if len(cand) >= 2:
        C = np.array([[G.cross_corr(a, b) for b in cand] for a in cand]); M = 0.5 * (C + C.T)
        ev = np.sort(np.linalg.eigvalsh(M))[::-1]
        rng = np.random.default_rng(0); nulls = []
        for _ in range(2000):
            N = rng.normal(0, 1 / math.sqrt(neff), (len(cand), len(cand))); nulls.append(np.linalg.eigvalsh(0.5 * (N + N.T)).max())
        thr = np.quantile(nulls, 1 - 0.00135)
        out.append(f'  eigen of sym C: {" ".join(f"{e:.3f}" for e in ev)}  null max (1-alpha) {thr:.3f} -> rank {(ev > thr).sum()}')
    if verbose: print('\n'.join(out))
    return out

if __name__ == '__main__':
    for d in sys.argv[1:]:
        analyse(d)
        print()

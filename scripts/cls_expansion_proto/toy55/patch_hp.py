#!/usr/bin/env python3
"""Adds the 'hp <A>' board option (shells coarser than hp excluded from the likelihood: infinite noise) to the single-process
klem board and to the distributed master and worker. usage: patch_hp.py SRC_DIR (src/main/commanders/simple)"""
import sys, os
d = sys.argv[1]
def patch(fn, reps):
    p = os.path.join(d, fn); s = open(p).read()
    if 'cfg%hp' in s or "case('hp')" in s: print('already patched', fn); return
    for o, n, c in reps:
        assert s.count(o) == c, (fn, o[:70], s.count(o)); s = s.replace(o, n)
    open(p, 'w').write(s); print('patched', fn)
cut_d = "if( cfg%hp > 0. .and. real(bd%shell(q)) < real(box) * smpd / cfg%hp ) sig2(q) = 1.e30"
patch('simple_klem_distr.f90', [
    ("pose_step = 3., pose_prior = 4., devw = 0.", "pose_step = 3., pose_prior = 4., devw = 0., hp = 0.", 1),
    ("                case('pnoise');        read(tok(2),*) cfg%pnoise\n", "                case('pnoise');        read(tok(2),*) cfg%pnoise\n                case('hp');            read(tok(2),*) cfg%hp\n", 1),
    ("\n                    if( bd%shell(q) <= kmax ) sig2(q) = real(sig2_sh(bd%shell(q)))\n",
     "\n                    if( bd%shell(q) <= kmax ) sig2(q) = real(sig2_sh(bd%shell(q)))\n                    " + cut_d + "\n", 1),
    ("\n                            if( bd%shell(q) <= kmax ) sig2(q) = real(sig2_sh(bd%shell(q)))\n",
     "\n                            if( bd%shell(q) <= kmax ) sig2(q) = real(sig2_sh(bd%shell(q)))\n                            " + cut_d + "\n", 1)])
patch('simple_commanders_klem_board.f90', [
    ("!!   specmatch 0|1 (every grid cell's radial amplitude spectrum matched to its half's consensus after the M-step)\n",
     "!!   specmatch 0|1 (every grid cell's radial amplitude spectrum matched to its half's consensus after the M-step)\n!!   hp 0 (A; > 0: shells coarser than hp carry no likelihood, so image background/offset noise cannot select particles)\n", 1),
    ("        real     :: pose_step, pose_prior, pgain, devw, fsc_s, pmed\n", "        real     :: pose_step, pose_prior, pgain, devw, fsc_s, pmed, hp\n", 1),
    ("pscale = 0; specmatch = 0;", "pscale = 0; specmatch = 0; hp = 0.;", 1),
    ("                case('pnoise');        read(tok(2),*) pnoise\n", "                case('pnoise');        read(tok(2),*) pnoise\n                case('hp');            read(tok(2),*) hp\n", 1),
    ("                    if( bd%shell(q) <= kmax ) sig2(q) = sig2_sh(bd%shell(q))\n",
     "                    if( bd%shell(q) <= kmax ) sig2(q) = sig2_sh(bd%shell(q))\n                    if( hp > 0. .and. real(bd%shell(q)) < real(box) * smpd / hp ) sig2(q) = 1.e30\n", 1)])

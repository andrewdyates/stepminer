#!/usr/bin/python
from __future__ import division
"""Compute all-pairs boolean implications.

Workflow script.

NOTE: if dual matrix comparison,
  rows of M1 map to rows of class matrix C
  rows of M2 map to cols of class matrix C

EXAMPLE:
/usr/bin/time python $HOME/pymod/stepminer/script.py m1_fname=/nfs/01/osu6683/gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/Methyl_correct_aligned.tab m2_fname=/nfs/01/osu6683/gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/mRNA_correct_aligned.sample.tab outdir=/nfs/01/osu6683
"""
import sys
from matrix_io import *
import script_thresholds
from __init__ import *
from lab_util import *

def main(m1_fname=None, m2_fname=None, percentile=.05, outdir="", cleanup=True, conf_t=2/3, stat_min=3, error_max=0.1):
  percentile, conf_t, stat_min, error_max = float(percentile), float(conf_t), float(stat_min), float(error_max)
  if outdir:
    make_dir(outdir)
  if isinstance(cleanup, str) and cleanup.lower() in ('f', 'false', 'none'):
    cleanup = False

  thresholds1 = script_thresholds.get_thresholds(fname=m1_fname, outdir=outdir, cleanup=cleanup)
  if m2_fname:
    thresholds2 = script_thresholds.get_thresholds(fname=m2_fname, outdir=outdir, cleanup=cleanup)

  D1 = load(m1_fname)
  c1 = get_bound(D1['M'], percentile=percentile)
  print "Bound for M1 based on 2*std of %.2f percentile std: %.4f" % (percentile*100, c1)
  if m2_fname:
    D2 = load(m2_fname)
    assert np.size(D1['M'],0) == np.size(D2['M'],0)
    c2 = get_bound(D2['M'], percentile=percentile)
    print "Bound for M2 based on 2*std of %.2f percentile std: %.4f" % (percentile*100, c2)

  # For all pairs, classify pair
  if not m2_fname:
    # Self compare.
    raise Exception, "Self comparison not yet implemented"  # TODO
  else:
    # Compare rows between two matrices.
    m,n = np.size(D1['M'],0), np.size(D2['M'],0)
    C = np.zeros((m,n), dtype=np.int)
    print "Created %d by %d classification matrix." % (m,n)
    for i, v1 in enumerate(D1['M']):
      for j, v2 in enumerate(D2['M']):
        c = classify_pair(*compress(v1, v2))
        if c:
          C[i,j] = CLASSES_I[c]

  # Save classifications
  out_fname = os.path.basename(m1_fname)
  if m2_fname:
    out_fname += ".%s" % m2_fname
  out_fname += ".stepminer.classes.pkl"
  out_fpath = os.path.join(outdir, out_fname)
  print "Saving stepminer boolean classes as %s..." % (out_fpath)
  save(C, out_fpath)
  print "Saved!"
  print "NOTE: classifications are enumerated integers as follows:"
  print CLASSES_I

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)

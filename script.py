#!/usr/bin/python
from __future__ import division
"""Compute all-pairs boolean implications.

Workflow script.

NOTE: if dual matrix comparison,
  rows of M1 map to rows of class matrix C
  rows of M2 map to cols of class matrix C

EXAMPLE:
/usr/bin/time python $HOME/pymod/stepminer/script.py m1_fname=/nfs/01/osu6683/gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/Methyl_correct_aligned.tab m2_fname=/nfs/01/osu6683/gse15745_nov2012_experiments/gse15745_aligned_matrices_nov2/mRNA_correct_aligned.tab outdir=/nfs/01/osu6683 cleanup=False
"""
import sys
from matrix_io import *
import script_thresholds
from __init__ import *
from lab_util import *

REPORT_ROW = 10

def main(m1_fname=None, m2_fname=None, percentile=.05, outdir="", cleanup=False, conf_t=2/3, stat_min=3, error_max=0.1):
  percentile, conf_t, stat_min, error_max = float(percentile), float(conf_t), float(stat_min), float(error_max)
  if outdir:
    make_dir(outdir)
  if isinstance(cleanup, str) and cleanup.lower() in ('f', 'false', 'none'):
    cleanup = False

  D1 = load(m1_fname)
  thresholds1 = script_thresholds.get_thresholds(fname=m1_fname, outdir=outdir, cleanup=cleanup, D=D1)
  if m2_fname:
    D2 = load(m2_fname)
    assert np.size(D1['M'],1) == np.size(D2['M'],1), "%s %s" % (D1['M'].shape,D2['M'].shape)
    thresholds2 = script_thresholds.get_thresholds(fname=m2_fname, outdir=outdir, cleanup=cleanup, D=D2)

  # ??????????????????????????????
  # HACK! Truncate D1 and D2
  # ??????????????????????????????
  print "TRUNCATING D1 and D2 for testing purposes!"
  print "original:", D1['M'].shape, D2['M'].shape
  D1['M'] = D1['M'][0:1000,:]
  D2['M'] = D2['M'][0:1000,:]
  print "new:", D1['M'].shape, D2['M'].shape

  c1 = get_bound(D1['M'], percentile=percentile)
  print "Bound for M1 based on 2*std of %.2f percentile std: c1=2*std=%.4f" % (percentile*100, c1)
  if m2_fname:
    c2 = get_bound(D2['M'], percentile=percentile)
    print "Bound for M2 based on 2*std of %.2f percentile std: c2=2*std=%.4f" % (percentile*100, c2)

  # For all pairs, classify pair
  if not m2_fname:
    # Self compare.
    raise Exception, "Self comparison not yet implemented"  # TODO
  else:
    C = classify_all_dual(D1, D2, thresholds1, thresholds2, c1, c2)

  # Save classifications
  out_fname = os.path.basename(m1_fname)
  if m2_fname:
    out_fname += ".%s" % os.path.basename(m2_fname)
  out_fname += ".stepminer.classes.pkl"
  out_fpath = os.path.join(outdir, out_fname)
  print "Saving stepminer boolean classes as %s..." % (out_fpath)
  save(C, out_fpath)
  print "Saved!"
  print "NOTE: classifications are enumerated integers as follows:"
  print CLASSES_I

  # Report stats for classification
  print "Class counts for %d items in (%d,%d) class matrix." % \
      (np.size(C), C.shape[0], C.shape[1])
  for name, x in CLASSES_I.items():
    cnt = np.sum(C==x)
    pct = cnt/np.size(C)*100
    print "%s (%d): %d (%.4f%%)" % (name, x, cnt, pct)
  print "Complete!"


if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  print args
  main(**args)

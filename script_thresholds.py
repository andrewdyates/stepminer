#!/usr/bin/python
"""Call Stepminer on gene expression .tab matrix to find midpoint thresholds.
Can be called as a script or imported as a sub routine (as in script.py).
Uses old Sahoo Java implementation.

EXAMPLE:
  python script_thresholds.py outdir=$HOME fname=example/mRNA10.tab
"""
from __init__ import *
import sys

def main(fname=None, outdir="", cleanup=False):
  """Return path to threshold per in-order row ID output file."""
  print "Loading %s..." % fname
  D = load(fname)
  row_thresholds = get_thresholds(fname=fname, outdir=outdir, cleanup=cleanup, D=D)
  thresh_out_fpath = os.path.join(outdir, get_thresh_out_fname(os.path.basename(fname)))
  print "Saving row thresholds in input matrix row order in file %s" % (thresh_out_fpath)
  fp = open(thresh_out_fpath, 'w')
  for rowid in D['row_ids']:
    fp.write("%s\t%f\n" % (rowid, row_thresholds[rowid]))
  fp.close()
  print "Saved %s." % (thresh_out_fpath)
  print "Script complete."
  return thresh_out_fpath, D

def get_thresholds(fname=None, outdir="", cleanup=False, overwrite=False, D=None):
  """Compute and return stepminer thresholds per row."""
  assert fname
  if not D:
    print "Loading matrix %s..." % fname
    D = load(fname)
  assert "row_ids" in D
  
  # check to see if stepminer input already exists
  expected_fname = os.path.join(outdir, get_stepminer_in_fname(os.path.basename(fname)))
  if not overwrite and os.path.exists(expected_fname):
    print "%s stepminer input exists and overwrite False, do not recreate file." % (expected_fname)
    Msort = np.sort(D['M'])
    stepmine_in_fname = expected_fname
  else:
    print "%s does not exist." % expected_fname
    print "Formatting .tab input matrix %s for stepminer input." % (fname)
    stepmine_in_fname, Msort = format_for_stepminer(fname, outdir, D=D)
    print "Loaded and sorted .tab matrix (%s by %s)" % D['M'].shape
    print "Created stepminer input file %s." % (stepmine_in_fname)

  expected_fname = get_stepminer_out_fname(stepmine_in_fname)
  if not overwrite and os.path.exists(expected_fname):
    print "%s stepminer output exists and overwrite is False, do not run stepminer and recreate file." % (expected_fname)
    stepmine_out_fname = expected_fname
  else:
    print "Running stepminer..."
    stepmine_out_fname = run_stepminer(stepmine_in_fname, verbose=True)
  print "Stepminer results saved in file %s." % (stepmine_out_fname)
  print "Parsing stepminer results..."
  row_thresholds = transform_output(open(stepmine_out_fname), D['row_ids'], Msort)
  print "Found %d thresholds from %s." % (len(row_thresholds), stepmine_out_fname)
  if cleanup:
    print "Deleting intermediate files..."
    os.remove(stepmine_in_fname); print "Deleted", stepmine_in_fname
    os.remove(stepmine_out_fname); print "Deleted", stepmine_out_fname
  return row_thresholds #, Msort
  

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  main(**args)

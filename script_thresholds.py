#!/usr/bin/python
"""Call stepminer on gene expression .tab matrix to find midpoint thresholds.

python script_thresholds.py outdir=$HOME/Desktop fname=example/mRNA10.tab
"""
from __init__ import *
import sys

def main(fname=None, outdir="", cleanup=False):
  """Return path to threshold per in-order row ID output file."""
  row_thresholds, row_ids = get_thresholds(fname=fname, outdir=outdir, cleanup=cleanup)
  thresh_out_fpath = os.path.join(outdir, get_thresh_out_fname(os.path.basename(fname)))
  print "Saving row thresholds in input matrix row order in file %s" % (thresh_out_fpath)
  fp = open(thresh_out_fpath, 'w')
  for rowid in row_ids:
    fp.write("%s\t%f\n" % (rowid, row_thresholds[rowid]))
  fp.close()
  print "Saved %s." % (thresh_out_fpath)
  print "Script complete."
  return thresh_out_fpath

def get_thresholds(fname=None, outdir="", cleanup=False):
  assert fname
  print "Formatting .tab input matrix %s for stepminer input." % (fname)
  stepmine_in_fname, D = format_for_stepminer(fname, outdir)
  print "Loaded and sorted .tab matrix (%s by %s)" % D['M'].shape
  print "Created stepminer input file %s." % (stepmine_in_fname)
  print "Running stepminer..."
  stepmine_out_fname = run_stepminer(stepmine_in_fname, verbose=True)
  print "Stepminer results saved in file %s." % (stepmine_out_fname)
  print "Parsing stepminer results..."
  row_thresholds = transform_output(open(stepmine_out_fname), D['row_ids'], D['M'])
  if cleanup:
    print "Deleting intermediate files..."
    os.remove(stepmine_in_fname); print "Deleted", stepmine_in_fname
    os.remove(stepmine_out_fname); print "Deleted", stepmine_out_fname
  return row_thresholds, D['row_ids']
  

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  main(**args)

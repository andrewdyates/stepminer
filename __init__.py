#!/usr/bin/python
"""
python script.py outdir=$HOME/Desktop fname=$HOME/Desktop/sahootest/mRNA10.tab
"""
from matrix_io import *
import subprocess


CMD = "java -jar %s -t OneStep %%s -o %%s" % (os.path.abspath("stepminer-1.1.jar"))

def stepminer_fname(fname, outdir):
  return os.path.join(outdir, os.path.basename(fname) + ".stepmine.pcl")

def stepminer_out_fname(fname, outdir):
  return os.path.join(outdir, os.path.basename(fname) + ".stepmine.pcl.ann")

def format_for_stepminer(fname, outdir):
  """Save a copy of an expression matrix for stepminer."""
  D = load(fname)
  for i in xrange(D['M'].shape[0]):
    D['M'][i] = sorted(D['M'][i])
  fp = open(stepminer_fname(fname, outdir), 'w')
  fp.write('YORF\tNAME\tGWEIGHT\t')
  fp.write("\t".join(map(str,xrange(D['M'].shape[1]))))
  fp.write('\nEWEIGHT\t\t\t')
  fp.write("\t".join(["1"]*D['M'].shape[1]))
  fp.write('\n')
  hacked_row_ids = ["%s\t%s\t1" % (s,s) for s in D['row_ids']]
  save(D['M'], fp, ftype='txt', row_ids=hacked_row_ids)
  fp.close()
  return stepminer_fname(fname, outdir)

def run_stepminer(fname, outdir):
  in_fname = format_for_stepminer(fname, outdir)
  out_fname = stepminer_out_fname(fname, outdir)
  cmd = CMD % (in_fname, out_fname)
  print cmd
  subprocess.call(cmd, shell=True)
  fp = open(stepminer_out_fname(fname, outdir))
  # TODO: extract rowid, p-value, threshold idx
  # TODO: save p-value, row_id, threshold idx, threshold value in original row order
  fp.next()
  for i, line in enumerate(fp):
    row = line.split('\t')
    rowid, pv, t = row[0], row[4], row[8]

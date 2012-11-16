#!/usr/bin/python
from __future__ import division
from matrix_io import *
import subprocess

CMD = "time java -jar %s -t OneStep %%s -o %%s" % (os.path.abspath("stepminer-1.1.jar"))

def get_stepminer_in_fname(fname):
  return fname + ".stepmine.pcl"

def get_stepminer_out_fname(fname):
  return fname + ".ann"

def get_thresh_out_fname(fname):
  return fname + ".stepmine_thresh.txt"

def format_for_stepminer(fname, outdir):
  """Save a copy of an expression matrix for stepminer, return loaded, sorted matrix."""
  D = load(fname)
  for i in xrange(D['M'].shape[0]):
    D['M'][i] = sorted(D['M'][i])
  stepminer_in_fpath = os.path.join(outdir, get_stepminer_in_fname(os.path.basename(fname)))
  fp = open(stepminer_in_fpath, 'w')
  fp.write('YORF\tNAME\tGWEIGHT\t')
  fp.write("\t".join(map(str,xrange(D['M'].shape[1]))))
  fp.write('\nEWEIGHT\t\t\t')
  fp.write("\t".join(["1"]*D['M'].shape[1]))
  fp.write('\n')
  hacked_row_ids = ["%s\t%s\t1" % (s,s) for s in D['row_ids']]
  save(D['M'], fp, ftype='txt', row_ids=hacked_row_ids)
  fp.close()
  return stepminer_in_fpath, D

def run_stepminer(in_fpath, verbose=True):
  out_fpath = get_stepminer_out_fname(in_fpath)
  cmd = CMD % (in_fpath, out_fpath)
  if verbose:
    print cmd
  subprocess.call(cmd, shell=True)
  return out_fpath

def transform_output(fp, row_ids, M):
  """Parse p-value, row_id, threshold idx, from .ann stepminer output file;
       Return threshold value in original row order.
  Args:
    fp: [*str] open file pointer to stepminer .ann results
    row_ids: [str] of row IDs in original matrix order
    M: np.array of per-row column-sorted original matrix; rows in original order
  Returns:
    obj?
  """
  fp.next() # skip first line
  rows = {}
  rowid_idx = dict([(s, i) for i, s in enumerate(row_ids)])
  for i, line in enumerate(fp):
    row = line.split('\t')
    rowid, t = row[0], row[8]
    # associate row ID with threshold value
    rows[rowid] = M[rowid_idx[rowid],int(t)]
  return rows
    
def get_bound(M, percentile=.05):
  """Return threshold +/- bound based on stdev percentile."""
  c = np.percentile(np.std(M, 1), percentile)
  return c

def classify_pair(v1, v2, t1, c1, t2, c2, conf_t=2/3):
  """Boolean classify pair of vectors.

  Args:
    v1, v2: [num] of vector pair to classify
    t1, t2: step threshold value for v1 and v2 respectively
    c1, c2: threshold confidence interval for v1 and v2 respectively
    conf_t: maximum ratio of values in threshold interval to classify
  Returns:
    ?
  """
  assert len(v1) == len(v2)
  low_conf1 = np.sum((v1 >= t1-c1) & (v1 <= t1+c1))
  low_conf2 = np.sum((v2 >= t2-c2) & (v2 <= t2+c2))
  if low_conf1/len(v1) > conf:
    return "Unclassifed"
  if low_conf2/len(v2) > conf:
    return "Unclassifed"
  
  a00 = np.sum(v1<t1-c1 & v2<t2-c2)
  a01 = np.sum(v1<t1-c1 & v2>t2+c2)
  a10 = np.sum(v1>t1+c1 & v2<t2-c2)
  a11 = np.sum(v1>t1+c1 & v2>t2+c2)
  total = a00 + a01 + a10 + a11
  assert low_conf1 + low_conf2 + total == len(v1) 

  expected = (a00+a01)*(a00+a10)/total
  observed = a00
  statistic = (expected-observed)/np.sqrt(expected)

  error_rate = 1/2 * (a00/(a00+a01) + a00/(a00+a10))

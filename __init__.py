#!/usr/bin/python
from __future__ import division
from matrix_io import *
import subprocess
import numpy as np

CMD = "time java -jar %s -t OneStep %%s -o %%s" % (os.path.abspath("stepminer-1.1.jar"))

CLASSES = {
  'UNL': set(), # Unspecified Nonlinear
  'HIH': set(['10']), # High Implies High
  'PC': set(['10', '01']), # Positive Correlation
  'LIL': set(['01']), # Low Implies Low
  'HIL': set(['11']), # High Imples Low
  'NC': set(['00', '11']), # Negative Correlation
  'LIH': set(['00']) # Low Implies High
}
CLASSES_I = dict((c, i+1) for i, c in enumerate(CLASSES))


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

def compress(v1, v2):
  """Return v1, v2 pair after removing any pairs with masked values."""
  assert len(v1)==len(v2)
  if hasattr(v1, 'mask'):
    m1 = v1.mask
  else:
    m1 = np.zeros(len(v1), dtype=np.bool)
  if hasattr(v2, 'mask'):
    m2 = v2.mask
  else:
    m2 = np.zeros(len(v2), dtype=np.bool)
  if not any(m1) and not any(m2):
    return v1, v2
  else:
    if not isinstance(v1, np.ndarray):
      v1v = np.array(v1)
    else:
      v1 = v1v
    if not isinstance(v2, np.ndarray):
      v2v = np.array(v2)
    else:
      v2v = v2
    return v1v[~(m1|m2)], v2v[~(m1|m2)]
  

def get_bound(M, percentile=.05):
  """Return threshold +/- bound based on stdev percentile."""
  c = np.percentile(np.std(M, 1), percentile)
  return c*2

def classify_pair(v1, v2, t1, c1, t2, c2, conf_t=2/3, stat_min=3, error_max=0.1):
  """Boolean classify pair of vectors.

  Args:
    v1, v2: [num] of vector pair to classify
    t1, t2: step threshold value for v1 and v2 respectively
    c1, c2: threshold confidence interval for v1 and v2 respectively
    conf_t: maximum ratio of values in threshold interval to classify
  Returns:
    str of CLASS in ('UNL', 'HIH', 'PC', 'LIL', 'HIL', 'NC', 'LIH')
  """
  assert len(v1) == len(v2)
  low_conf1 = np.sum((v1 >= t1-c1) & (v1 <= t1+c1))
  low_conf2 = np.sum((v2 >= t2-c2) & (v2 <= t2+c2))
  if low_conf1/len(v1) > conf_t:
    return "Unclassifed"
  if low_conf2/len(v2) > conf_t:
    return "Unclassifed"

  counts = {
    '00': np.sum(v1<t1-c1 & v2<t2-c2),
    '01': np.sum(v1<t1-c1 & v2>t2+c2),
    '10': np.sum(v1>t1+c1 & v2<t2-c2),
    '11': np.sum(v1>t1+c1 & v2>t2+c2),
    }
  sparse = {}
  for q in counts:
    stat, err = test_sparse(q,counts)
    sparse[q] = (stat >= stat_min and err <= error_max)
  # Classify based on which quads are sparse.
  all_quads = set(sparse.keys())
  for cls, quads in CLASSES.items():
    if all(sparse[q] for q in quads) and \
          all(not sparse[q] for q in all_quads-quads):
      return cls
  # This should never return None; if so, there is some bug.
  return None
  
def test_sparse(q, D):
  """Return if quadrant q is sparse given quad density counts dict D"""
  assert q in ('00', '01', '10', '11')
  assert set(D.keys()) == set(('00', '01', '10', '11'))
  total = sum(D.keys())
  
  expected = (D[q] + D[q[0]+toggle(q[1])]) * (D[q] + D[toggle(q[0])+q[1]]) / total
  observed = D[q]
  statistic = (expected - observed) / np.sqrt(expected)
  error = 0.5 * (D[q] / (D[q] + D[q[0]+toggle(q[1])]) + D[q] / (D[q] + D[toggle(q[0])+q[1]]))
  return statistic, error

def toggle(c):
  if c == '0': return '1'
  elif c == '1': return '0'
  else: return None

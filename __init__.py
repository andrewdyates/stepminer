#!/usr/bin/python
from __future__ import division
from matrix_io import *
import subprocess
import numpy as np

this_dir = os.path.dirname(os.path.abspath(__file__))
CMD = "/usr/bin/time java -jar %s -t OneStep %%s -o %%s" % (os.path.join(this_dir, "stepminer-1.1.jar"))

# Notation:
#  quadrant 'yx' corresponds to (row_var, col_var) where (0,1)=>(low,high)
#  e.g., '01' means "row var Y low, col var X high"
#  Class name always noted in terms of "X implies Y" or "Col_var implies Row_var"
# --------------------
# Class Symbol: strict set of sparse quadrants
CLASSES = {
  'UNL': set(),            # Unspecified Nonlinear: special case
  'HIH': set(['01']),      # High Implies High: (Low Y, High X) sparse
  'PC': set(['01', '10']), # Positive Correlation: (Low Y, High X) and (High Y, Low X) sparse
  'LIL': set(['10']),      # Low Implies Low: (Low X, Low Y) sparse
  'HIL': set(['11']),      # High Imples Low: (High Y, High X) sparse
  'NC': set(['00', '11']), # Negative Correlation: (Low Y, Low X) and (High Y, High X) sparse
  'LIH': set(['00'])       # Low Implies High: (Low Y, Low X) sparse
}
CLASSES_I = {
  'UNL': 0,
  'HIH': 1,
  'PC': 2,
  'LIL': 3,
  'HIL': 4,
  'NC': 5,
  'LIH': 6,
}
ALL_QUADS = set(('00', '01', '10', '11'))


def set_quantile(D,c,thresholds):
  """Return matrix of (-1, 0, 1) of element row quantile."""
  Mq = np.zeros(D['M'].shape, dtype=np.int8)
  for i, row in enumerate(D['M']):
    t = thresholds[D['row_ids'][i]]
    low = D['M'][i]<(t-c)
    high = D['M'][i]<(t+c)
    # low = 0 (default), high = 1, nonsig = 4 (not high or low)
    Mq[i,high] = 1 
    Mq[i,~high&~low] = 4
  return Mq

def test_sparse(q, D):
  """Return if quadrant q is sparse given quad density counts dict D
  q[0] is Y quantile (rows in matrix)
  q[1] is X quantile (cols in matrix)
  e.g., q<-'01' means "quadrant of low Y, high X"
  
  Args:
    q: str of quadrant in (yx) = ('00', '01', '10', '11')
    D: {str:int} of element counts in each quadrant
  """
  assert q in ALL_QUADS
  assert set(D.keys()) == ALL_QUADS

  Yq = D[q] + D[q[0]+toggle(q[1])] 
  Xq = D[q] + D[toggle(q[0])+q[1]]
  total = sum(D.values())
  expected = Yq * Xq / total
  observed = D[q]
  # --- Joint probability statistic.
  statistic = (expected - observed) / np.sqrt(expected)
  # --- Error Rate
  error = 0.5 * (observed/Yq + observed/Xq)
  return statistic, error

def toggle(c):
  if c == '0': return '1'
  else: return '0'

def classify_from_sparse(sparse):
  """Return classification from "quadrant is sparse" declarations."""
  assert set(sparse.keys()) == ALL_QUADS
  for cls, quads in CLASSES.items():
    if all(sparse[q] for q in quads) and \
          all(not sparse[q] for q in ALL_QUADS-quads):
      return cls
  # This should never return None; if so, there is some bug.
  return None
  
def classify_all_dual(D1, D2, thresholds1, thresholds2, c1, c2, conf=2/3, stat_min=3, error_max=0.1):
  """Boolean classify all rows in two matrix objects.

  Args:
    D1: matrix dict of row variables (in scatterplot: y axis)
    D2: matrix dict of col variables (in scatterplot: x axis)
    thresholds1: {str=>float} of per row midpoint thresholds for D1, indexed by row name
    thresholds2: {str=>float} of per row midpoint thresholds for D2, indexed by row name
    c1: float of confidence interval, c1
    c2: float of confidence interval, c2
    conf: fraction of points in confidence interval to make bool class insignificant
    stat_min: num of min "z-score" like test for quad sparse significance (~ pv<=0.001)
    error_max: num of maximum error for quad sparse significance
  """
  # low = 0 (default), high = 1, nonsig = 4 (not high or low)
  m1,m2,n = np.size(D1['M'],0), np.size(D2['M'],0), np.size(D1['M'],1)
  Mq1 = set_quantile(D1,c1,thresholds1)
  Mq2 = set_quantile(D2,c2,thresholds2)

  # Q: are enough elements in either high or low regions?
  Q1 = np.sum(Mq1 == 4,1) <= n*nonsig
  Q2 = np.sum(Mq2 == 4,1) <= n*nonsig
  # Create classification matrix. By default, all entries are 0 => UNL
  C = np.zeros((m1,m2), dtype=np.int8)
  
  # For each pair of rows such that both have enough points in either high and low
  #   i: row index (y axis); j: col index (x axis)
  Mq1=Mq1*2 # shift all bits one to the left in Mq1 (Y)
  for i in np.nonzero(Q1):
    for j in np.nonzero(Q2):
      y, x = Mq1[i], Mq2[j]
      r = y+x # (0,1)->(high,low); msb: y, lsb: x
      counts = {
       '00': np.sum(r==0),
       '01': np.sum(r==1),
       '10': np.sum(r==2),
       '11': np.sum(r==3)
       }
      sparse = {}
      for q in counts:
        stat, err = test_sparse(q,counts)
        sparse[q] = (stat >= stat_min and err <= error_max)
      cls = classify_from_sparse(sparse)
      if cls is None:
        print "WARNING: Could not classify row %d=>%s, col %d=>%s (idx from 0)" % \
            (i, D1['row_ids'][i], j, D2['row_ids'][j])
      else:
        C[i,j] = CLASSES_I[cls]
  return C

def get_stepminer_in_fname(fname):
  return fname + ".stepmine.pcl"

def get_stepminer_out_fname(fname):
  return fname + ".ann"

def get_thresh_out_fname(fname):
  return fname + ".stepmine_thresh.txt"

def format_for_stepminer(fname, outdir, D=None):
  """Save a copy of an expression matrix for stepminer, return loaded, sorted matrix.

  Returns:
    str of filepath to stepminer-formated input
    np.array of per row sorted array
  """
  if not D:
    D = load(fname)
  Msort = np.sort(D['M'])
  stepminer_in_fpath = os.path.join(outdir, get_stepminer_in_fname(os.path.basename(fname)))
  fp = open(stepminer_in_fpath, 'w')
  fp.write('YORF\tNAME\tGWEIGHT\t')
  fp.write("\t".join(map(str,xrange(D['M'].shape[1]))))
  fp.write('\nEWEIGHT\t\t\t')
  fp.write("\t".join(["1"]*D['M'].shape[1]))
  fp.write('\n')
  hacked_row_ids = ["%s\t%s\t1" % (s,s) for s in D['row_ids']]
  # Save sorted copy of matrix to fp with headed already added.
  save(Msort, fp, ftype='txt', row_ids=hacked_row_ids)
  fp.close()
  return stepminer_in_fpath, Msort

def run_stepminer(in_fpath, verbose=True):
  out_fpath = get_stepminer_out_fname(in_fpath)
  cmd = CMD % (in_fpath, out_fpath)
  if verbose:
    print cmd
  subprocess.call(cmd, shell=True)
  return out_fpath

def transform_output(fp, row_ids, Msort):
  """Parse p-value, row_id, threshold idx, from .ann stepminer output file;
       Return threshold value in original row order.
  Args:
    fp: [*str] open file pointer to stepminer .ann results
    row_ids: [str] of row IDs in original matrix order
    M: np.array of per-row column-sorted original matrix; rows and cols in original order
  Returns:
    {str: num} of rowID=>threshold
  """
  fp.next() # skip first line
  rows = {}
  rowid_idx = dict([(s, i) for i, s in enumerate(row_ids)])
  for i, line in enumerate(fp):
    row = line.split('\t')
    rowid, t = row[0], int(row[8])-3 # account for first three meta rows
    # associate row ID with threshold value
    idx = rowid_idx[rowid]
    try:
      rows[rowid] = (Msort[idx,t] + Msort[idx,t+1]) / 2.0
    except IndexError:
      print i, rowid, t, rowid_idx[rowid]
      raise
  return rows

def get_bound(M, percentile=.05):
  """Return threshold +/- bound based on stdev percentile."""
  c = np.percentile(np.std(M, 1), percentile)
  return c*2


# DEPRECIATED
# ========================================
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
  if not isinstance(v1, np.ndarray):
    v1 = np.array(v1)
  if not isinstance(v2, np.ndarray):
    v2 = np.array(v2)
  low_conf1 = np.sum((v1 >= t1-c1) & (v1 <= t1+c1))
  low_conf2 = np.sum((v2 >= t2-c2) & (v2 <= t2+c2))
  if low_conf1/len(v1) > conf_t:
    return "UNL"
  if low_conf2/len(v2) > conf_t:
    return "UNL"

  counts = {
    '00': np.sum((v1<t1-c1) & (v2<t2-c2)),
    '01': np.sum((v1<t1-c1) & (v2>t2+c2)),
    '10': np.sum((v1>t1+c1) & (v2<t2-c2)),
    '11': np.sum((v1>t1+c1) & (v2>t2+c2)),
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
  

def get_quantile_thresholds(D,c,thresholds):
  """Return get highest-low, lowest-high indices per row."""
  Msort = np.sort(D['M'])
  L = np.zeros(Msort.shape[0],2, dtype=np.int)
  for i, row in enumerate(Msort):
    t = thresholds[D['row_ids'][i]]
    L[i] = np.searchsorted(row, (t-c,t+c))
  return L

#!/usr/bin/python
# not quite... needs to be in same format
from matrix_io import *
FNAME = "mRNA10.tab"
D = load(FNAME)
for i in xrange(D['M'].shape[0]):
  D['M'][i] = sorted(D['M'][i])

fp = open("mRNA10.sort.pcl", 'w')
fp.write('YORF\tNAME\tGWEIGHT\t')
fp.write("\t".join(map(str,xrange(D['M'].shape[1]))))
fp.write('\nEWEIGHT\t\t\t')
fp.write("\t".join(["1"]*D['M'].shape[1]))
fp.write('\n')
hacked_row_ids = ["%s\t%s\t1" % (s,s) for s in D['row_ids']]
save(D['M'], fp, ftype='txt', row_ids=hacked_row_ids)
fp.close()

#! /bin/env python

# Author: Elliot Drabek
# will report # of contigs, cumulative length, and n50 of an assembly

from sys import *
from drabek import read_fasta

###############################################################################

lengths = []
total_length = 0
for id, seq in read_fasta():
  lengths.append(len(seq))
  total_length += len(seq)

half_total = total_length / 2.0

lengths.sort()
total_so_far = 0
num_so_far = 0
for length in reversed(lengths):
  total_so_far += length
  num_so_far += 1
  if total_so_far >= half_total:
    n50 = length
    break

#print >> stderr, 'contigs\tbases\tn50'
print '\t'.join(map(str, (len(lengths), total_length, n50)))

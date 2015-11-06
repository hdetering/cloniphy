#!/usr/bin/env python3
import dendropy
import dendropy.calculate.popgenstat as pgs
import sys

if len(sys.argv) < 2:
    print("usage: %s alignment.fasta" % sys.argv[0], file=sys.stderr)
    sys.exit(0)

# read alignment
alignment_fasta = sys.argv[1]
char_mat = dendropy.DnaCharacterMatrix.get(file=open(alignment_fasta), schema='fasta')

# calculate stats
n_seg_sites = pgs.num_segregating_sites(char_mat)

# report
print("segregating sites:\t%d" % n_seg_sites)

#!/usr/bin/env python3
from __future__ import division, print_function
import os
import sys
from Bio import SeqIO

if __name__ == '__main__':
    fn_seq = sys.argv[1]
    if not os.path.exists(fn_seq):
        print("[ERROR] File does not exist: ".format(fn_seq), file=sys.stderr)

    print("[INFO] Scanning genome file: {}".format(fn_seq), file=sys.stderr)

    count = {}
    glen = 0
    for rec in SeqIO.parse(fn_seq, 'fasta'):
        print("[INFO] Scanning chromosome '{}'".format(rec.id), file=sys.stderr)

        pos = 0
        trinuc = ''
        for n in str(rec.seq):
            # we don't count 'N's
            if n == 'N':
                trinuc = ''
                continue
            # make sure three nucleotides have been consumed yet
            if len(trinuc) == 3:
                if trinuc not in count:
                    count[trinuc] = 1
                else:
                    count[trinuc] += 1
                trinuc = "{}{}".format(trinuc[1:], n)
            else:
                trinuc = "{}{}".format(trinuc, n)
            pos += 1
        glen += pos

    # print trinucleotide frequencies
    for trinuc, c in sorted(count.items()):
        print('{0},{1:.8f}'.format(trinuc, c/glen))

#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert mutation matrix (MM) to SNV-to-cluster mapping CNV.
#
# Background: 
#   In the MM it is not clear, which mutations are SNVs and which CNVs.
#   Mutations are only identifiable by their column index which is not ideal
#   for lookup.
# Usage:
#   awk -f mm_to_snv.awk somatic.vcf clones.nex
# Input:
#   - somatic.vcf: contains SNV ids
#   - clones.nex: clone tree with clusters and MM
#     The column index in the MM corresponds to SNV id in somatic.vcf)
# Output:
#   id_snv,chrom,pos,cluster
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-27
#------------------------------------------------------------------------------

# print CSV header
BEGIN {
  printf("id_mut,chrom,pos,id_cluster\n");
}
# read SNV ids
NR == FNR && !/^#/ {
  n_snvs++;
  id_snv = $3
  chrom[id_snv] = $1;
  pos[id_snv] = $2;
  next;
}
# turn on mutation matrix parsing
/Matrix$/ {
  mm = 1;
  next;
}
# skip comments (enclosed by quare brackets)
mm && /\[.*\]/ {
  next;
}
# turn off mutation  matrix parsing
mm && (NF != 2) {
  mm = 0;
}
# parse mutation matrix
mm && (NF == 2) {
  for (i=1; i<=length($2); i++) {
    if (substr($2, i, 1) == "1")
      if ("s"i in chrom)
        printf("s%d,%s,%s,%s\n", i, chrom["s"i], pos["s"i], $1);
  }
}
END {
  #printf("parsed %d snvs\n", length(chrom));
}

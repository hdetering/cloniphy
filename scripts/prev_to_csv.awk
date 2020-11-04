#!/usr/bin/awk
# vim: syntax=awk tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert prevalence matrix (samples x clones) to long CSV.
#
# Output:
#   id_cluster,id_sample,freq
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-10-11
#------------------------------------------------------------------------------

BEGIN {
  FS  = ",";
  OFS = ",";
  print "id_cluster", "id_sample", "freq";
}
# get column names from header row
NR == 1 {
  for (i=2; i<=NF; i++) {
    cname[i] = $i;
  }
  next;
}
# parse data rows
{
  for (i=2; i<=NF; i++) {
    print cname[i], $1, $i;
  }
}

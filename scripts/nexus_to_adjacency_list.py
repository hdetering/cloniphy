#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert tree in NEXUS notation to adjacency list.
#
# Usage:
#   python3 nexus_to_adjacency_list.py tree.nex
# Output:
#   from,to
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-27
#------------------------------------------------------------------------------

import argparse
from nexus import NexusReader
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Convert NEXUS tree to adjacency list (CSV).')
  parser.add_argument('tree', type=argparse.FileType('r'), help='Path to NEXUS tree file.')
  
  args = parser.parse_args()
  return args

def parse_tree(filename):
  # load tree from NEXUS file
  n = NexusReader.from_file(filename)
  nwk = n.trees[0].split('=')[1].strip()
  # convert string representation to Tree (format=1: allow for internal node names)
  t = Tree(nwk, format=1)
  adj = []
  for node in t.traverse('preorder'):
    for child in node.children:
      adj.append((node.name, child.name))
  
  return adj

if __name__ == '__main__':
  args = parse_args()
  adj = parse_tree(args.tree.name)
  for node, child in adj:
    print('{},{}'.format(node, child))
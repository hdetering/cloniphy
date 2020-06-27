#!/usr/bin/env python3
# vim: syntax=python tabstop=2 shiftwidth=2 expandtab
# coding: utf-8
#------------------------------------------------------------------------------
# Convert tree in Newick notation to adjacency list.
#
# Usage:
#   awk -f tree_to_adjacency_list.awk tree.nwk
# Output:
#   from,to
#------------------------------------------------------------------------------
# author   : Harald Detering
# email    : harald.detering@gmail.com
# modified : 2020-06-27
#------------------------------------------------------------------------------

import argparse
from ete3 import Tree

def parse_args():
  parser = argparse.ArgumentParser(description='Convert Newick tree to adjacency list (CSV).')
  parser.add_argument('tree', type=argparse.FileType('r'), help='Path to Newick tree file.')
  
  args = parser.parse_args()
  return args

def parse_tree(filename):
  # load tree from file (format=1: allow for internal node names)
  t = Tree(filename, format=1)
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

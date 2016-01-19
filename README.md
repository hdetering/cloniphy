# What is CloniPhy?

CloniPhy is a simulator to generate genomic sequences based on a clone phylogeny.

# How does the simulator work?

The starting point is a reference genome plus (optionally) a set of germline mutations that are used as a personal genome. A genealogy of clones can be given as a Newick tree or generated randomly (given the number of clones). A set of mutations are introduced to the clone tree and subsequently applied to the personal genome, yielding an individual genome for each clone.

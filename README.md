# What is CloniPhy?

CloniPhy is a simulator that generates genomic sequencing data based on a clone phylogeny.

# How does the simulator work?

1. The starting point is a *reference genome* plus (optionally) a set of *germline mutations* that is used to generate a personal genome. 
2. A *genealogy of clones* can be given as a Newick tree or generated randomly (given the number of clones). 
3. A set of *somatic mutations* are introduced to the clone tree and subsequently applied to the personal genome, yielding an individual genome for each clone.
4. *Sequencing reads* are sampled from each clone genome.

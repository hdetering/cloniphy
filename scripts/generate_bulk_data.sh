#!/bin/bash
set -x

# parameters for simulation
N_CLONES=5     # number of clones
FREQ=(0.5 0.25 0.15 0.05 0.05) # relative clone frequencies
N_MUT=100      # number of mutations
N_TRANS=20     # number of transforming mutations
PURITY=0.7     # cellularity (proportion of tumor cells in the bulk)
COVERAGE=30    # target sequencing coverage
READ_LEN=150   # read length
#N_READS=100000 # total number of reads
REF="/home/harry/code/cloniphy/data/genome_1Mb.fa" # healthy genome

# simulate clone genomes
bin/cloniphy -c $N_CLONES -f ${FREQ[@]} -m $N_MUT -t $N_TRANS -r $REF

# set up read simulator
ART="/home/harry/src/art_bin_ChocolateCherryCake/art_illumina"
ART_PARAMS="-sam -p -l $READ_LEN -ss HS25 -m 500 -s 10"
# simulate sequencing reads (in proportion to clone frequency)
$ART $ART_PARAMS -i $REF -f $(bc -l <<< "(1-$PURITY)*$COVERAGE") -o healthy_reads > art.log
for i in $(seq 1 $N_CLONES); do
  my_idx=$(($i-1))
  my_cvg=$(bc -l <<< "$PURITY * ${FREQ[$my_idx]} * $COVERAGE / 2")
  my_pfx=$(printf "clone%02d" $i)
  $ART $ART_PARAMS -i ${my_pfx}_genome.fa -f $my_cvg -o ${my_pfx}_reads >> art.log
done

# pool reads artificially
cat *_reads1.fq > pool_reads1.fq
cat *_reads2.fq > pool_reads2.fq

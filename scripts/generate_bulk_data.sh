#!/bin/bash
set -x

# parameters for simulation
N_CLONES=5     # number of clones
FREQ=(0.5 0.25 0.15 0.05 0.05) # relative clone frequencies
N_MUT=2000     # number of mutations
N_TRANS=100    # number of transforming mutations
PURITY=0.7     # cellularity (proportion of tumor cells in the bulk)
COVERAGE=50    # target sequencing coverage
READ_LEN=150   # read length
#N_READS=100000 # total number of reads
REF="/home/harry/code/cloniphy/data/hs37d5_chr20.fa" # healthy genome
REF_VCF="/home/harry/code/cloniphy/data/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.HG00128.genotypes.vcf"

# simulate clone genomes
date
echo "Simulating clone genomes..."
bin/cloniphy -c $N_CLONES -f ${FREQ[@]} -m $N_MUT -t $N_TRANS -r $REF -v $REF_VCF > clonify.log 2>&1
date

# set up read simulator
ART="/home/harry/src/art_bin_ChocolateCherryCake/art_illumina"
ART_PARAMS="-sam -na -p -l $READ_LEN -ss HS25 -m 500 -s 10"
# simulate sequencing reads (in proportion to clone frequency)
echo "Simulating sequencing reads..."
$ART $ART_PARAMS -i "healthy_genome.fa" -f $(bc -l <<< "(1-$PURITY)*$COVERAGE / 2") -o healthy_reads > art.log
for i in $(seq 1 $N_CLONES); do
  my_idx=$(($i-1))
  my_cvg=$(bc -l <<< "$PURITY * ${FREQ[$my_idx]} * $COVERAGE / 2")
  my_pfx=$(printf "clone%02d" $i)
  $ART $ART_PARAMS -i ${my_pfx}_genome.fa -f $my_cvg -o ${my_pfx}_reads >> art.log
done
date

# compress output
echo "Compressing output files"
gzip healthy_reads?.fq clone??_reads.fq
for f in *.sam; do
  samtools -bS $f > ${f%sam}bam
done
date

# pool reads artificially
zcat healthy_reads1.fq clone??_reads1.fq | gzip pool_reads1.fq.gz
zcat healthy_reads2.fq clone??_reads2.fq | gzip pool_reads2.fq.gz

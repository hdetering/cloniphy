#!/bin/bash

# paths to input data
REF="/home/harry/code/cloniphy/data/input/hs37d5.fa"
READS1="/home/harry/code/cloniphy/data/output/pool_reads1.fq.gz"
READS2="/home/harry/code/cloniphy/data/output/pool_reads2.fq.gz"
REF_INDELS="/home/harry/code/cloniphy/data/input/1000G_phase1.indels.b37.vcf"
# tools to be used
GATK="java -jar /home/harry/src/gatk/GenomeAnalysisTK.jar"
PICARD="java -jar /home/harry/src/picard-tools-1.141/picard.jar"

# index reference
if [ ! -f $REF.bwt ]; then
  bwa index $REF
fi
# map reads to reference
if [ ! -e pool.bam ]; then
  bwa mem -M -R "@RG\tID:1\tSM:HG00128\tPL:Illumina\tLB:Sim1\tPU:0" \
  $REF $READS1 $READS2 \
  | samtools view -bS - > pool.bam
fi
# sort bam
if [ ! -e pool.s.bam ]; then
  samtools sort -O bam pool.bam > pool.s.bam
fi

# mark duplicates
if [ ! -e pool.dedup.bam ]; then
  $PICARD MarkDuplicates \
    INPUT=pool.s.bam \
    OUTPUT=pool.dedup.bam \
    METRICS_FILE=pool.dedup.metrics \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT
fi

# prepare reference for use with GATK
if [ ! -f ${REF%fa}dict ]; then
  $PICARD CreateSequenceDictionary R=$REF O=${REF%fa}dict
fi
if [ ! -f $REF.fai ]; then
  samtools faidx $REF
fi

# realign around indels
if [ ! -e target.intervals.list ]; then
  $GATK -T RealignerTargetCreator -R $REF \
    -o target.intervals.list \
    -known $REF_INDELS
fi
if [ ! -e pool.realigned.bam ]; then
  $GATK -T IndelRealigner -R $REF \
    -I pool.dedup.bam \
    -targetIntervals target.intervals.list \
    -o pool.realigned.bam \
    -known $REF_INDELS
fi

# filter reads by mapping quality
samtools view -b -q 40 pool.realigned.bam > pool.MQfiltered.bam

# Adjust quality scores
#$GATK -T BaseRecalibrator \
#    -R $REF \
#    -I pool.MQfiltered.bam \
#    -o pool.MQFiltered.table \
#    -knownSites DbSNP_138.vcf
$GATK -T PrintReads \
    -R $REF \
    -I pool.MQfiltered.bam \
    -o pool.recalibrated.bam \
    -BQSR pool.MQfiltered.table

# call variants
$GATK -T UnifiedGenotyper -R $REF \
  -I pool.recalibrated.bam \
  -stand_call_conf 30.0 -stand_emit_conf 10.0 \
  -glm BOTH \
  -dcov 2500 -mbq 20 \
  -o pool.raw.vcf

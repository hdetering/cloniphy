#!/bin/bash
set -e
# prefix must be provided as parameter (e.g. "sample1")
PFX=$1

# paths to input data
REF="/home/harry/code/cloniphy/data/multisample/ref.fa"
READS="/home/harry/code/cloniphy/data/multisample/$PFX.fq"
# tools to be used
GATK="java -jar /home/harry/src/gatk/GenomeAnalysisTK.jar"
PICARD="java -jar /home/harry/src/picard-tools-2.5.0/picard.jar"

# index reference
if [ ! -f $REF.bwt ]; then
  bwa index $REF
fi
# map reads to reference
#if [ ! -e ${PFX}.m.bam ]; then
#  bwa mem -M -R "@RG\tID:${PFX}\tSM:${PFX}\tPL:Illumina\tLB:Sim1\tPU:0" \
#  $REF $READS \
#  | samtools view -bS - > ${PFX}.m.bam
#fi
# compress simulated SAM file
if [ ! -e ${PFX}.bam ]; then
  samtools view -bS ${PFX}.sam > ${PFX}.bam
fi

# sort bam
if [ ! -e ${PFX}.s.bam ]; then
  samtools sort -O bam ${PFX}.bam > ${PFX}.s.bam
fi

# mark duplicates
if [ ! -e ${PFX}.dedup.bam ]; then
  $PICARD MarkDuplicates \
    INPUT=${PFX}.s.bam \
    OUTPUT=${PFX}.dedup.bam \
    METRICS_FILE=${PFX}.dedup.metrics \
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

# filter reads by mapping quality
if [ ! -f ${PFX}.MQfilt.bam ]; then
  samtools view -b -q 40 ${PFX}.dedup.bam > ${PFX}.MQfilt.bam
  samtools index ${PFX}.MQfilt.bam
  $GATK -T DepthOfCoverage -R ref.fa -o ${PFX}.cvg -I ${PFX}.MQfilt.bam -pt readgroup
fi

# call variants
$GATK -T HaplotypeCaller -R $REF \
  -I ${PFX}.MQfilt.bam \
  -stand_call_conf 30.0 -stand_emit_conf 10.0 \
  -mbq 20 \
  -o ${PFX}.vcf

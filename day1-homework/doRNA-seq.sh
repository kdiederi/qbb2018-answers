#!/bin/bash

GENOME=~/qbb2018-answers/genomes/BDGP6
ANNOTATION=~/data/genomes/BDGP6.Ensembl.81.gtf

for SAMPLE in SRR072893 SRR072903 SRR072905 SRR072915
do
    echo "making directory"
    mkdir $SAMPLE
    cd $SAMPLE
    echo "running fastqc"
    fastqc ~/data/rawdata/${SAMPLE}.fastq
    echo "running hisat2"
    hisat2 -x $GENOME -U ~/data/rawdata/${SAMPLE}.fastq -S ${SAMPLE}_map.sam
    echo "running samtools sort"
    samtools sort -o ${SAMPLE}_map.bam ${SAMPLE}_map.sam
    echo "running samtools index"
    samtools index -b ${SAMPLE}_map.bam
    echo "running stringtie"
    stringtie ${SAMPLE}_map.bam -p -e -G $ANNOTATION -B -o ${SAMPLE}_map.gtf
    cd ../
done

#!/bin/bash

set -e


##BAM and indexing##
#Creation of sorted BAM file with filtered alignments for yeast
#Determination of continuous intervals for filtered alignments for references
#Creation of sorted BAM file with mapped-unmapped alignments for yeast
head -n 2 results/sam/reads_t.sam > results/sam/transposon_intersect_2.sam
cat results/sam/transposon_intersect.sam >> results/sam/transposon_intersect_2.sam

head -n 17 results/sam/reads.sam > results/sam/yeast_intersect_2.sam
cat results/sam/yeast_intersect.sam >> results/sam/yeast_intersect_2.sam

samtools view -S -b results/sam/transposon_intersect_2.sam > results/bam/transposon_intersect_2.bam
samtools view -S -b results/sam/yeast_intersect_2.sam > results/bam/yeast_intersect_2.bam

bedtools bamtobed -i results/bam/transposon_intersect_2.bam > results/bed/transposon_intersect_2.bed
bedtools bamtobed -i results/bam/yeast_intersect_2.bam > results/bed/yeast_intersect_2.bed

bedtools merge -i results/bed/transposon_intersect_2.bed > docs/transposon_intervals_2.txt
bedtools merge -i results/bed/yeast_intersect_2.bed > docs/yeast_intervals_2.txt

samtools sort results/bam/yeast_intersect_2.bam -o results/bam/yeast_intersect_2_sorted.bam
samtools index results/bam/yeast_intersect_2_sorted.bam results/bam/yeast_intersect_2_sorted.bam.bai

samtools sort results/bam/combined.bam -o results/bam/combined_sorted.bam
samtools index results/bam/combined_sorted.bam results/bam/combined_sorted.bam.bai


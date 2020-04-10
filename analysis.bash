#!/bin/bash

set -e

##folder structure##
#In a base directory, the folders "data", "docs", "results" and "src" are created
#In the "results" folder, the folders "bam", "bed", "sam" and "fastqc" are created
#The compressed sequencing data is kept in the "data" folder
#This code file is kept in the "src" folder
#Ensure you are in the base directory before running this bash file using "bash src/analysis.bash"


##quality check##
#Unzips the sequence data file and runs fastqc on each read file
#Results are stored in "results/fastqc" folder and summary is stored in "docs" folder
tar -xvzf data/A0166851X.tgz
fastqc data/*.fq -o results/fastqc

for filename in results/fastqc/*.zip; do
       unzip -d results/fastqc $filename
done

cat results/fastqc/*/summary.txt > docs/fastqc_summary.txt


##mapped-unmapped alignment##
#Builds index for reference sequences and filters out mate-pairs where...
#one read is mapped to the reference but the other cannot be
#Determination of intervals that the aligned reads span with results stored in "docs" folder
#SAM files containing the filtered mate-pairs are created for each reference
bowtie2-build data/sacCer3.fa data/yeast_genome

export BOWTIE2_INDEXES=$(pwd)/data

bowtie2 -x yeast_genome -p 4 -1 data/reads_1.fq -2 data/reads_2.fq\
               -S results/sam/reads.sam

samtools view -S -b results/sam/reads.sam > results/bam/reads.bam
samtools sort results/bam/reads.bam -o results/bam/reads_sorted.bam
samtools view -b -F 4 -f 8 results/bam/reads_sorted.bam > results/bam/mapped_unmapped.bam
samtools view -b -F 8 -f 4 results/bam/reads_sorted.bam > results/bam/unmapped_mapped.bam
samtools merge -f results/bam/combined.bam results/bam/mapped_unmapped.bam results/bam/unmapped_mapped.bam
bedtools bamtobed -i results/bam/combined.bam > results/bed/combined.bed
bedtools merge -i results/bed/combined.bed > docs/yeast_intervals.txt

bowtie2-build data/ty5_6p.fa data/transposon

export BOWTIE2_INDEXES=$(pwd)/data

bowtie2 -x transposon -p 4 -1 data/reads_1.fq -2 data/reads_2.fq\
               -S results/sam/reads_t.sam

samtools view -S -b results/sam/reads_t.sam > results/bam/reads_t.bam
samtools sort results/bam/reads_t.bam -o results/bam/reads_sorted_t.bam
samtools view -b -F 4 -f 8 results/bam/reads_sorted_t.bam > results/bam/mapped_unmapped_t.bam
samtools view -b -F 8 -f 4 results/bam/reads_sorted_t.bam > results/bam/unmapped_mapped_t.bam
samtools merge -f results/bam/combined_t.bam results/bam/mapped_unmapped_t.bam results/bam/unmapped_mapped_t.bam
bedtools bamtobed -i results/bam/combined_t.bam > results/bed/combined_t.bed
bedtools merge -i results/bed/combined_t.bed > docs/transposon_intervals.txt

samtools view results/bam/combined.bam > results/sam/combined.sam
samtools view results/bam/combined_t.bam > results/sam/combined_t.sam

echo "Do you wish to continue?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) bash src/analysis2.bash; break;;
        No ) exit;;
    esac
done
echo "completed"


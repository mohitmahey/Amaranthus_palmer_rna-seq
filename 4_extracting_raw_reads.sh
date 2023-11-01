#!/bin/bash
# this script is used to extract raw reads from the aligned files

# setting the working directory
data="/mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum"

# loading required modules
module load GCC/11.3.0
module load SAMtools/1.16.1

cd ${data}

list=$(cat samples.list)

for file in $list;
do
  samtools idxstats -@ 64 ${data}/sorted_aligned_transcripts/${file}.bam > ${data}/raw_counts/${file}_counts.txt
  samtools stats -@ 64 ${data}/sorted_aligned_transcripts/${file}.bam > ${data}/transcripts_stats/${file}_stats

done







#!/bin/bash
# this script will take the cleaned reads and align it to genome
# using hisat2 aligner, which is a true alinger

# setting the working directory
data="/mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum"

# laoding required modules
module load GCC/11.3.0
module load SAMtools/1.16.1
module load hisat2/2.1.0

# building the genome index 
hisat2-build -p 64 ${data}/genome/Amaranthus_palmeri_reference.transcripts.fa \
${data}/genome/Amaranthus_palmeri

cd ${data}

list=$(cat samples.list)
for file in $list;
do
# using hisat2 aligner to align paired-end reads to genome, and producing SAM file
  	hisat2 -p 64 --quiet \
               -x ${data}/genome/Amaranthus_palmeri \
               -1 ${data}/cleaned_files/${file}_R1.fastq \
               -2 ${data}/cleaned_files/${file}_R2.fastq \
               -S ${data}/aligned_transcripts/${file}.sam 

# using smatools to sort and convert sam to bam, as sam file are very large
        samtools view -hb -@ 64 ${data}/aligned_transcripts/${file}.sam  | samtools sort -@ 64 \
        -o ${data}/sorted_aligned_transcripts/${file}.bam
 
# making index of the aligned files
        samtools index ${data}/sorted_aligned_transcripts/${file}.bam

# finally removing the sam file to save disk space
        rm ${data}/aligned_transcripts/${file}.sam
done


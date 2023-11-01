#!/bin/bash
# this script will run the BWA to align 

# setting the working dir
data="/mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum"

# changing to wd
cd ${data}

# loading the required modules
module load GCC/11.3.0
module load SAMtools
module load GCC/5.4.0-2.26  OpenMPI/1.10.3
module load BWA/0.7.17
module load picard/2.25.0-Java-11

# getting the samples name in a variable
list=$(cat ${data}/samples.list) 

# first step is to index the trancriptome
bwa index ${data}/genome/Amaranthus_palmeri_reference.transcripts.fa

# to align the files we will use bwa mem script
for file in ${list};
do
  bwa mem -t 64 -P -R "@RG\tID:${file}\tSM:${file}\tPL:ILLUMINA\tLB:${file}\tPU:RTSF" -o ${data}/bwa_transcripts/${file}_aligned.sam \
   ${data}/genome/Amaranthus_palmeri_reference.transcripts.fa \
   ${data}/cleaned_files/${file}_R1.fastq ${data}/cleaned_files/${file}_R2.fastq 

# converting sam to bam and sorting them
 samtools view -hb -@ 64 ${data}/bwa_transcripts/${file}_aligned.sam | samtools sort -@ 64  \
         -o ${data}/bwa_aligned/${file}_sorted.bam
samtools index ${data}/bwa_aligned/${file}_sorted.bam

# removing the sam file
rm ${data}/bwa_transcripts/${file}_aligned.sam
done

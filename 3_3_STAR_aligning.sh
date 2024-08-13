#!/bin/bash
# this script will take the cleaned reads and align it to genome
# using hisat2 aligner, which is a true alinger

# setting the working directory
data="/mnt/gs21/scratch/maheymoh/ksu_amrathus_rish"

# laoding required modules
module load GCC/11.3.0
module load STAR/2.7.10b

# building the genome index 
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ${data}/genome \
--genomeFastaFiles ${data}/genome/Amaranthus_palmeri_reference.transcripts.fa 

cd ${data}

#making directories
mkdir ${data}/star_alignments

list=$(cat samples.list)
for file in $list;
do
# runnning star aligner
STAR --runThreadN 16 --genomeDir ${data}/genome --readFilesIn ${data}/cleaned_files/${file}_R1.fastq ${data}/cleaned_files/${file}_R2.fastq \
 --outFileNamePrefix ${data}/star_alignments/${file}

# using smatools to sort and convert sam to bam, as sam file are very large
  #      samtools view -hb -@ 64 ${data}/star_alignments/${file}.sam  | samtools sort -@ 64 \
   #     -o ${data}/sorted_aligned_transcripts/${file}.bam
 
# making index of the aligned files
    
#  samtools index ${data}/sorted_aligned_transcripts/${file}.bam
# finally removing the sam file to save disk space
 #       rm ${data}/aligned_transcripts/${file}.sam
done


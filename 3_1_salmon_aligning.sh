#!/bin/bash
# this script aligns the transcripts using salmon

# setting the working directory
data="/mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum"

# loading required modules
module load GCC/10.2.0  OpenMPI/4.0.5
module load Salmon/1.8.0

cd ${data}

# building the index of transcritome to be used by salm
salmon index -t ${data}/genome/Amaranthus_palmeri_reference.transcripts.fa \
 -i ${data}/genome/AmaPa_index -k 31

cd ${data}
list=$(cat samples.list) 

for file in $list;
do
# running the salmon alignment
salmon --no-version-check \
       quant -i ${data}/genome/AmaPa_index -l A \
             -1 ${data}/cleaned_files/${file}_R1.fastq \
             -2 ${data}/cleaned_files/${file}_R2.fastq \
             --output ${data}/salmon_transcripts/${file}_quant.sf \
             -p 64
done









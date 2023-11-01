#!/bin/bash

# this script is first step of trimming and qc of raw data
# fastp is used for this 


# setting the working directory
data="/mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum"

# setting the conda env
source /mnt/home/maheymoh/miniconda3/bin/activate

# this has fastp
conda activate variants

# changing directory to raw_data files
cd ${data}/raw_data

# making directory for cleaned reads
mkdir ${data}/cleaned_files
mkdir ${data}/cleaned_files/qc


# this is to make list of the raw files
printf '%s\n' * | grep '_R1\_001\.fastq\.gz' > samples.list

# this just keep the first name and removes everything else
sed -i 's/\(.*\)_R1\_001\.fastq\.gz\1//' samples.list

list=$(cat samples.list)

for file in ${list};
do
   fastp -i ${file}_R1_001.fastq.gz  -I ${file}_R2_001.fastq.gz -o ${data}/cleaned_files/${file}_R1.fastq -O ${data}/cleaned_files/${file}_R2.fastq -w 64 \
   --dedup --failed_out ${data}/cleaned_files/fail.fq -j ${data}/cleaned_files/qc/${file}_fastp.json -h ${data}/cleaned_files/qc/${file}_fastp.html 
done

mv ${data}/raw_data/samples.list ${data}/samples.list




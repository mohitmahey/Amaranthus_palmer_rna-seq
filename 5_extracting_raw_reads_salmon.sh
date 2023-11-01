#!/bin/bash

# this script is to get raw reads from the salmon 
# output folders and put it in a .txt file.

# wd
data="/mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum/salmon_transcripts"

# changing to the salmon_transcripts folder
cd ${data}

# each folder is named based on the list.sample file
# each folder has file "quant.sf" with name and raw reads
# name is first column and raw reads is 4th column

# getting the sample list
list=$(cat /mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum/samples.list)

# make a for loop with name of list as input

for file in $list;
do
# cd into folders
  cd ${data}/${file}_quant.sf

# use awk '{print $1, $3}' file.txt to get columns
   awk '{print $1, $5}' quant.sf > ${data}/transcripts_reads_count/${file}_raw_counts.txt
done








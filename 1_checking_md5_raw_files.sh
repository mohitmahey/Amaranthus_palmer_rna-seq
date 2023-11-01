#!/bin/bash

# script to check the integrity of the raw files
# md5 sequences is provided by the company and we will be 
# running our own md5 to match the integrity

# we have file and then .md5 extension
# EASY, just use "md5sum --check" with command line tool

# change dir to the raw data folder
cd /mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum/raw_data
ls -1 *.md5 | while read line; do cat $line ; done | md5sum -c > md5_integrity_check.txt

cd /mnt/gs21/scratch/maheymoh/amranthus_KSU_dr_jugulum/raw_data_2
ls -1 *.md5 | while read line; do cat $line ; done | md5sum -c > md5_integrity_check_2.txt

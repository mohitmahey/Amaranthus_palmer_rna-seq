# Amaranthus_palmer_rna-seq

Differential gene expression anlaysis and Weighted co-expression analysis of Amarathus palmer. 

# Experimental Design

Susceptible and Resistant Amaranthus palmer was cloned and treated with five different hericides - Atrazine, 24D, chlorosulfuron. mesotrione.
The RNA was extracted and was sequenced using illumina hi-seq 150 bp sequencing.

# Algining the raw RNA-seq reads to the transcriptome

The raw reads was aligned to the transcriptome of amranthis palmer, to get the read counts data that was used further for the DGE and WGCNA analysis

# Tools used

1.fastp - for removing adapter sequences and quality check of raw fastq files Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

2.Hisat2 - for aligning the raw reads to the genome (Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019)

3.SAMTools - for converting and handling SAM/BAM files (Heng Li and others, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352)

4.EdgeR - for making comparisons between different populations and getting differentially expressed genes. http://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

# Procedure 

## 1 - Checking the integrity of the raw fastq files
The raw fastq file as received from the company, comes along with a hash file that can be used for checking integrity of files and to confirm that we received the files withour any loss of
data. For this we used md5sum hash algorithm to confirm integrity.(1_checking_md5_raw_files.sh)

```
ls -1 *.md5 | while read line; do cat $line ; done | md5sum -c > md5_integrity_check.txt
```
here we are making a list of all the .md5 files and combining it into a singlelist that can be used to to perform md5sum check.

## 2 - Trimming adapters, cleaing bad reads and quality check

fastp is a command-line tool that can trim the adapter sequences, filter low quality reads and provide quality control files for the raw reads. it provides the quality control files in html and json format,
to see the quality of the sequences. (2_fastp_trimming.sh)

```
fastp -i ${file}_R1_001.fastq.gz  -I ${file}_R2_001.fastq.gz -o ${data}/cleaned_files/${file}_R1.fastq -O ${data}/cleaned_files/${file}_R2.fastq -w 64 \
   --dedup --failed_out ${data}/cleaned_files/fail.fq -j ${data}/cleaned_files/qc/${file}_fastp.json -h ${data}/cleaned_files/qc/${file}_fastp.html
```
-i = input of forward strand
-I = input of reverse strand
-o = name of forward strand output file, after filtering
-O = name of reverse strand output file, after filtering
-w = number of cores to be used (max it can use = 16)
--dedup = drop duplicated sequences
--failed_out = file where all the low-quality reads are stored
-j = file name of quality control file in json format
-h = file name of quality control file in html format

## 3 - aliging the raw reads to the transcriptome

After initial filtering, trimming, and quality control of the raw reads, the next step is to align the raw reads. We used HISAT2 to align the raw reads to the reference transcriptome.
for redundancy we also used bwa and salmon to see if it caused any significant changes to the final count data.





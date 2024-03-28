# Amaranthus_palmeri_rna-seq

Differential gene expression analysis and Weighted co-expression analysis of Amaranthus palmeri. 

# Experimental Design

Susceptible and Resistant Amaranthus palmeri was cloned and treated with five different herbicides - Atrazine, 2,4,-D, chlorosulfuron, mesotrione and lactofen.
The RNA was extracted and sequenced using Illumina hi-seq 150 bp sequencing.

# Aligning the raw RNA-seq reads to the transcriptome

The raw reads were aligned to the transcriptome of amaranthus palmeri, to get the read count data that was used further for the DGE analysis and WGCNA. 

# Tools used

1. fastp - for removing adapter sequences and quality check of raw fastq files Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560

2. Hisat2 - for aligning the raw reads to the genome (Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019)

3. SAMTools - for converting and handling SAM/BAM files (Heng Li and others, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352)

4. EdgeR - for making comparisons between different populations and getting differentially expressed genes. http://bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

5. WGCNA - for weighed co-expression anaylsis among different populations, Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559

   
# Procedure 

## 1 - Checking the integrity of the raw fastq files(1_checking_md5_raw_files.sh)
The raw fastq file as received from the company, comes along with a hash file that can be used for checking integrity of files and to confirm that we received the files without any loss of data. For this we used md5sum hash algorithm to confirm integrity.

```
ls -1 *.md5 | while read line; do cat $line ; done | md5sum -c > md5_integrity_check.txt
```
here we are making a list of all the .md5 files and combining it into a singlelist that can be used to to perform md5sum check.

## 2 - Trimming adapters, cleaing bad reads and quality check(2_fastp_trimming.sh)

fastp is a command-line tool that can trim the adapter sequences, filter low quality reads and provide quality control files for the raw reads. it provides the quality control files in html and json format, to see the quality of the sequences. 

```
fastp -i ${file}_R1_001.fastq.gz  -I ${file}_R2_001.fastq.gz -o ${data}/cleaned_files/${file}_R1.fastq -O ${data}/cleaned_files/${file}_R2.fastq -w 64 \
   --dedup --failed_out ${data}/cleaned_files/fail.fq -j ${data}/cleaned_files/qc/${file}_fastp.json -h ${data}/cleaned_files/qc/${file}_fastp.html
```
-i = input of forward strand \
-I = input of reverse strand \
-o = name of forward strand output file, after filtering \
-O = name of reverse strand output file, after filtering \
-w = number of cores to be used (max it can use = 16) \
--dedup = drop duplicated sequences \
--failed_out = file where all the low-quality reads are stored \
-j = file name of quality control file in json format \
-h = file name of quality control file in html format

## 3 - aligning the raw reads to the transcriptome (3_hisat2_align_sam_to_bam.sh)

After initial filtering, trimming, and quality control of the raw reads, the next step is to align the raw reads. We used HISAT2 to align the raw reads to the reference transcriptome.
for redundancy we also used bwa and salmon to see if it caused any significant changes to the final count data. 

To run the aligner, we need 
a. cleaned forward reads (made in previous using fastp) \ 
b. cleaned reverse reads (made in previous using fastp) \
c. Indexed reference transcriptome to be aligned. \

```
# building the genome index 
hisat2-build -p 64 ${data}/genome/Amaranthus_palmeri_reference.transcripts.fa \
${data}/genome/Amaranthus_palmeri

hisat2 -p 64 --quiet \
               -x ${data}/genome/Amaranthus_palmeri \
               -1 ${data}/cleaned_files/${file}_R1.fastq \
               -2 ${data}/cleaned_files/${file}_R2.fastq \
               -S ${data}/aligned_transcripts/${file}.sam
```

where 
-p = number of cores to use \
--quiet = to prevent printing to terminal, except sequences or serious errors \
-x = input indexed reference transcriptome \
-1 = forward reads \
-2 = reverse reads \
-S = output as SAM format (default output is stdout) \

### Sorting and aligning sorted SAM outputs

Hisat2 ouputs as SAM files. Those can be manipulated ussing SAMTools. For downstream analysis, we need to sort the SAM files and index them.
SAM were converted to BAM only to save disk space as SAM files can be very large.
```
# using smatools to sort and convert sam to bam, as sam file are very large
        samtools view -hb -@ 64 ${data}/aligned_transcripts/${file}.sam  | samtools sort -@ 64 \
        -o ${data}/sorted_aligned_transcripts/${file}.bam
 
# making index of the aligned files
        samtools index ${data}/sorted_aligned_transcripts/${file}.bam

# finally removing the sam file to save disk space
        rm ${data}/aligned_transcripts/${file}.sam
```

## Extracting raw reads from aligned files(4_extracting_raw_reads.sh)

To get the mapped reads from the aligned BAM files we use SAMtools idxstat
```
samtools idxstats -@ 64 ${data}/sorted_aligned_transcripts/${file}.bam > ${data}/raw_counts/${file}_counts.txt
```
this gives you a file with transcripts name, sequence length, mapped reads number and unmapped reads number
for more information, see [samtools idxtstat](https://www.htslib.org/doc/samtools-idxstats.html)

you can parse it out using awk to get only transcripts name and mapped reads number only as thats only you need for the DGE analysis

```
awk '{print $1, $3}' counts_${file}.txt > counts_reads_${file}.txt
```

After you have the counts file in the txt form, you can combine all of them in excel to get combine_counts file, that can 
be used as an input in the edgeR package to get Differentially expressed genes. Save it as .txt file that can be loaded in R later on.

# Diffential gene expression analysis using EdgeR

EdgeR was used to find differentially expressed genes.

## Dissecting the code

Firstly, the environment is cleaned and the required packages are loaded 

```
# clearing envs
rm(list=ls())

# installing the package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

# loading the package
library(edgeR)
```

The combined_counts data is loaded in to the R environment

```
# loading in the combined counts data 
x <- read.table("combined_counts_amrPa_ksu.txt", header = TRUE, row.names = "GeneID", sep = '\t')
```

For EdgeR, grouping of the treatments is needed. in this case as there are 5 treaments and a control, there will be 6 groups with 3 replicated each.
This grouping is important as to make comparisons later on.

```
# This is when comparing each with each, it follows the above classification
group <- factor(c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10))
```

then the groups are assigned to different treatments

```
# assigning groups to each group
y <- DGEList(counts = x,group = group)
y$samples
```

Filtering the low counts

```
# filtering out low reads from the reads counts
keep <- filterByExpr(y)
summary(keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples
```

Normalisation of the counts data to prevent bias based on the library size.

```
# Normalising the reads
y <- calcNormFactors(y)
y$samples
```

making  model design 
```
# making the design to do the comparisons and contrast, 10 groups are formed based on 
# grouping described in lines 13- 22
design <- model.matrix(~0+group, data = y$samples)
```
Running the model to find the differential expression

```
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
```

doing the comparisons , here the comparison between resistant population non-treated with susceptible population non-treated
```
qlf.R_NTvS_NT<- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0,0,0,0))
write.table(qlf.R_NTvS_NT, "R_NTvS_NT.txt")
```
























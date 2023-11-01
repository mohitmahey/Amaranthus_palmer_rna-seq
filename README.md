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

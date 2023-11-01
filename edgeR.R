# this code is to do DGEs of AmaPa counts from KSU
# clearing envs
rm(list=ls())

# installing the package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

# loading the package
library(edgeR)

# loading in the combined counts data 
x <- read.table("combined_counts_amrPa_ksu.txt", header = TRUE, row.names = "GeneID", sep = '\t')
View(x)

#1 = resistant NT, non-treated resistant samples
#2 = resistant ATR
#3 = resistant CHL
#4 = resistant MESO
#5 = resistant 24D
#6 = susceptible NT, non-treated susceptible samples
#7 = susceptible ATR
#8 = susceptible CHL
#9 = susceptible MESO
#10 = susceptible 24D

# This is when comparing each with each, it follows the above classification
group <- factor(c(1,1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10))

# This is when comparing SvR
# group <- factor(c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))

# This is when we are comparing ALl Sus vs each R
# group <- factor(c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))

# assigning groups to each group
y <- DGEList(counts = x,group = group)
y$samples

# filtering out low reads from the reads counts
keep <- filterByExpr(y)
summary(keep)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y$samples

# normalising the reads
y <- calcNormFactors(y)
y$samples

# making the design to do the comparisons and contrast, 10 groups are formed based on 
# grouping described in lines 13- 22
design <- model.matrix(~0+group, data = y$samples)

# design <- model.matrix(~0+group, data = y$samples) is same as 
# design <- model.matrix(~group-1, data = y$samples)

design

# writing the design table as txt file
write.table(design, "design.txt")


y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)

# sample clustering
plotMDS(y)

# doing the comparisons ------
# R vs S, non- treated meaning KCTRNT vs KSSNT, combing all the reps and NT and NT_2
# DGE that are constitutively expressed
qlf.R_NTvS_NT<- glmQLFTest(fit, contrast = c(1,0,0,0,0,-1,0,0,0,0))
write.table(qlf.R_NTvS_NT, "R_NTvS_NT.txt")
R_NTvS_NT <- read.table("R_NTvS_NT.txt")
View(R_NTvS_NT)

#----------ATR-------------------------
# R,T vs S,T, treated with ATR, meaning KCTRATR vs KSSATR
# DGE within bio types
qlf.R_ATRvS_ATR <- glmQLFTest(fit, contrast = c(0,1,0,0,0,0,-1,0,0,0))
write.table(qlf.R_ATRvS_ATR, "R_ATRvS_ATR.txt")
R_ATRvS_ATR <- read.table("R_ATRvS_ATR.txt")
View(R_ATRvS_ATR)

# R, NT vs R, T - treated with ATR, meaning KCTRNT vs KCTRATR
# DGE induced by herbicide
qlf.R_NTvR_ATR <- glmQLFTest(fit, contrast = c(1,-1,0,0,0,0,0,0,0,0))
write.table(qlf.R_NTvR_ATR, "R_NTvR_ATR.txt")
R_NTvR_ATR <- read.table("R_NTvR_ATR.txt")
View(R_NTvR_ATR)

# S, NT vs S, T - treated with ATR, meaning KSSNT vs KCTRATR
# DGE induced by herbicide in susceptible, but does not give resistance
qlf.S_NTvS_ATR <- glmQLFTest(fit, contrast = c(0,0,0,0,0,1,-1,0,0,0))
write.table(qlf.S_NTvS_ATR, "S_NTvS_ATR.txt")
S_NTvS_ATR <- read.table("S_NTvS_ATR.txt")
View(S_NTvS_ATR)

#----------CHL-------------------------
# R,T vs S,T, treated with CHL, meaning KCTRCHL vs KSSCHL
# DGE within bio types
qlf.R_CHLvS_CHL <- glmQLFTest(fit, contrast = c(0,0,1,0,0,0,0,-1,0,0))
write.table(qlf.R_CHLvS_CHL, "R_CHLvS_CHL.txt")
R_CHLvS_CHL <- read.table("R_CHLvS_CHL.txt")
View(R_CHLvS_CHL)

# R, NT vs R, T - treated with ATR, meaning KCTRNT vs KCTRATR
# DGE induced by herbicide
qlf.R_NTvR_CHL <- glmQLFTest(fit, contrast = c(1,0,-1,0,0,0,0,0,0,0))
write.table(qlf.R_NTvR_CHL, "R_NTvR_CHL.txt")
R_NTvR_CHL <- read.table("R_NTvR_CHL.txt")
View(R_NTvR_CHL)

# S, NT vs S, T - treated with ATR, meaning KSSNT vs KCTRATR
# DGE induced by herbicide in susceptible, but does not give resistance
qlf.S_NTvS_CHL <- glmQLFTest(fit, contrast = c(0,0,0,0,0,1,0,-1,0,0))
write.table(qlf.S_NTvS_CHL, "S_NTvS_CHL.txt")
S_NTvS_CHL <- read.table("S_NTvS_CHL.txt")
View(S_NTvS_CHL)

#----------Meso-------------------------
# R,T vs S,T, treated with ATR, meaning KCTRATR vs KSSATR
# DGE within bio types
qlf.R_MESOvS_MESO <- glmQLFTest(fit, contrast = c(0,0,0,1,0,0,0,0,-1,0))
write.table(qlf.R_MESOvS_MESO, "R_MESOvS_MESO.txt")
R_MESOvS_MESO <- read.table("R_MESOvS_MESO.txt")
View(R_MESOvS_MESO)

# R, NT vs R, T - treated with ATR, meaning KCTRNT vs KCTRATR
# DGE induced by herbicide
qlf.R_NTvR_MESO <- glmQLFTest(fit, contrast = c(1,0,0,-1,0,0,0,0,0,0))
write.table(qlf.R_NTvR_MESO, "R_NTvR_MESO.txt")
R_NTvR_MESO <- read.table("R_NTvR_MESO.txt")
View(R_NTvR_MESO)

# S, NT vs S, T - treated with ATR, meaning KSSNT vs KCTRATR
# DGE induced by herbicide in susceptible, but does not give resistance
qlf.S_NTvS_MESO <- glmQLFTest(fit, contrast = c(0,0,0,0,0,1,0,0,-1,0))
write.table(qlf.S_NTvS_MESO, "S_NTvS_MESO.txt")
S_NTvS_MESO <- read.table("S_NTvS_MESO.txt")
View(S_NTvS_MESO)

#----------24D-------------------------
# R,T vs S,T, treated with ATR, meaning KCTRATR vs KSSATR
# DGE within bio types
qlf.R_24DvS_24D <- glmQLFTest(fit, contrast = c(0,0,0,0,1,0,0,0,0,-1))
write.table(qlf.R_24DvS_24D, "R_24DvS_24D.txt")
R_24DvS_24D <- read.table("R_24DvS_24D.txt")
View(R_24DvS_24D)

# R, NT vs R, T - treated with ATR, meaning KCTRNT vs KCTRATR
# DGE induced by herbicide
qlf.R_NTvR_24D <- glmQLFTest(fit, contrast = c(1,0,0,0,-1,0,0,0,0,0))
write.table(qlf.R_NTvR_24D, "R_NTvR_24D.txt")
R_NTvR_24D <- read.table("R_NTvR_24D.txt")
View(R_NTvR_24D)

# S, NT vs S, T - treated with ATR, meaning KSSNT vs KCTRATR
# DGE induced by herbicide in susceptible, but does not give resistance
qlf.S_NTvS_24D <- glmQLFTest(fit, contrast = c(0,0,0,0,0,1,0,0,0,-1))
write.table(qlf.S_NTvS_24D, "S_NTvS_24D.txt")
S_NTvS_24D <- read.table("S_NTvS_24D.txt")
View(S_NTvS_24D)


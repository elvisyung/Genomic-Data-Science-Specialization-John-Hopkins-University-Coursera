# Week 4 - Module 4 Quiz (Questions: Rsamtools, oligo, limma, minfi, edgeR and DESEQ2 functions)

# Question 1
library(yeastRNASeq)
fastqFilePath <- system.file("reads", "wt_1_f.fastq.gz", package = "yeastRNASeq")

library(ShortRead)
library(Biostrings)

# read fastq file
reads <- readFastq(fastqFilePath)
DNAStringSet <- sread(reads)

# fraction of reads has an A in the 5th base
cm <- consensusMatrix(DNAStringSet, as.prob=TRUE, baseOnly=TRUE)
cm['A', 5]

# Output
# 0.363841

# Question 2
mean(as(quality(reads), "matrix")[,5])

# Output
# 28.93346

# Question 3
library(leeBamViews)
bamFilePath <- system.file("bam", "isowt5_13e.bam", package="leeBamViews")

library(Rsamtools)

bamFile <- BamFile(bamFilePath)

# focus on Scchr13, interval from 800,000 to 801,000
gr <- GRanges(seqnames = "Scchr13", ranges = IRanges(start = c(800000), end = c(801000)))
params <- ScanBamParam(which = gr, what = scanBamWhat())
aln <- scanBam(bamFile, param = params)

# find duplicates
sum(table(aln[[1]]$pos)) - sum(table(aln[[1]]$pos) == 1)

# Output
# 129

# Question 4
bpaths <- list.files(system.file("bam", package="leeBamViews"), pattern = "bam$", full=TRUE)

# focus on the novel transcribed regions
bamView <- BamViews(bpaths)
gr_nt <- GRanges(seqnames="Scchr13", ranges=IRanges(start = c(807762), end = c(808068)))
bamRanges(bamView) <- gr_nt
aln_nt <- scanBam(bamView)

# get sequences for each sample
alns <- lapply(aln_nt, function(xx) xx[[1]]$seq)

# calculate the average number of reads across 8 the samples
alns_len_sum = 0
for (i in 1:length(alns)){
  alns_len_sum = alns_len_sum + length(alns[i][[1]])
}
alns_len_sum / length(alns)

# Output
# 90.25

# Question 5
library(oligo)
library(GEOquery)

# get data
getGEOSuppFiles("GSE38792")
untar("GSE38792/GSE38792_RAW.tar", exdir = "GSE38792/CEL")

# read data
celfiles <- list.files("GSE38792/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)

# parse pData
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename)
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)), "OSA", "Control")

# find "8149273" probeset
normData <- rma(rawData)
loc <- match("8149273", rownames(normData))

mean(exprs(normData[loc,])[1:8])

# Output
# 7.02183

# Question 6
library(limma)

# use limma to fit between control group and OSA group
normData$group <- factor(normData$group)
design <- model.matrix(~normData$group)
fit <- lmFit(normData, design)
fit <- eBayes(fit)

# absolute value of logFC which has lowest P.value
abs(topTable(fit)$logFC[1])

# Output
# 0.7126484

# Question 7
fit_toptable <- topTable(fit)
de <- subset(fit_toptable, adj.P.Val < 0.05)
de

# Output
# <0 rows> (or 0-length row.names)

# Question 8 
library(minfi)
library(minfiData)

# get OpenSea loci in RGsetEx with preprocess
rgSet <- preprocessFunnorm(RGsetEx)
rg_opensea <- rgSet[c(getIslandStatus(rgSet) == "OpenSea")]

# get Beta value in both group
rg_beta <- getBeta(rg_opensea)
normal <- mean(rg_beta[, c(1,2,5)])
cancer <- mean(rg_beta[, c(3,4,6)])

# mean difference between normal and cancer group
normal - cancer

# Output
# 0.08863657

# Question 9
library(AnnotationHub)

# get Caco2 data
ah <- AnnotationHub()
ah <- subset(ah, species=="Homo sapiens")
ah_Caco2 <- query(ah, c("Caco2", "AWG"))
ah_Caco2 <- ah_Caco2[["AH22442"]]

CpG_450K <- granges(rgSet)

unique(findOverlaps(CpG_450K, ah_Caco2, type="within"))

# Outout
# queryLength: 485512 / subjectLength: 123048

# Question 10
library(DESeq2)
library(zebrafishRNASeq)

# get and parse data
data("zfGenes")
zf <- zfGenes[grep("^ERCC", rownames(zfGenes), invert = T), ]
zf <- as.matrix(zf)
colData <- DataFrame(sampleID = colnames(zf), group = as.factor(c("control", "control", "control", "treatment", "treatment", "treatment")))

# perform DESeq2
dds <- DESeqDataSetFromMatrix(zf, colData, design = ~ group)
dds <- DESeq(dds)

# find differentially expressed features
res <- results(dds)
sigRes <- subset(res, padj <= 0.05)
dim(sigRes)[1]

# Output
# 116


# Week 3 - Module 3 Quiz (Questions: ExpressionSet, GEOquery and biomaRT usuage and functions)

# Question 1
library(ALL)

data(ALL)
means(exprs(ALL[,5]))

# output 
# 5.629627

# Question 2
library(biomaRt)
library("hgu95av2.db")

# list Marts
mart <- useMart(host='feb2014.archive.ensembl.org', biomart = "ENSEMBL_MART_ENSEMBL")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)

# annotate each feature
feature_name <- featureNames(ALL)
annotation_ALL <- getBM(attributes=c("ensembl_gene_id","affy_hg_u95av2"), filters="affy_hg_u95av2", values=feature_name, mart=ensembl)

sum(table(annotation_ALL[,2])>1)

# Output 
# 1045

# Question 3
# list Attributes
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

# annotate autosomes
chrom <- c(1:22)
annotation_ALL_chr <- getBM(attributes=c("ensembl_gene_id", "affy_hg_u95av2", "chromosome_name"), filters=c("affy_hg_u95av2","chromosome_name"), values=list(feature_name, chrom), mart=ensembl)

sum(table(table(annotation_ALL_chr[,2])))

# Output
# 11016

# Question 4
library(minfiData)
library(minfi)

mean(getMeth(MsetEx)[,2])

# Output
# 7228.277

# Question 5
library(GEOquery)

eList <- getGEO("GSE788")
eData <- eList[[1]]

mean(exprs(eData)[,2])

# Output
# 756.432

# Question 6
library(airway)
library(GenomicRanges)

data(airway)
mean(airway$avgLength)

# Output
# 113.75

# Question 7
sum(assay(airway)[,3]>=1)

# Output
# 25699

# Question 8
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# exon data of txdb
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb_exons <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)

# transcripts on the autosome
autosome <- paste0("chr", c(1:22))
txdb_exons_autosome <- keepSeqlevels(txdb_exons, autosome, pruning.mode = "coarse")

# rename in NCBI format
txdb_ncbi <- mapSeqlevels(seqlevels(txdb_exons), "NCBI")
txdb_exons_ncbi <- renameSeqlevels(txdb_exons_autosome, txdb_ncbi)

dim(subsetByOverlaps(airway, txdb_exons_ncbi))[1]

# Output
# 26276

# Question 9 
sample_SRR1039508 <- airway[, 1]
sample_SRR1039508_autosome <- subsetByOverlaps(sample_SRR1039508, txdb_exons_ncbi)

autosome_reads <- sum(assay(sample_SRR1039508_autosome, "counts"))
total_reads <- sum(assay(sample_SRR1039508, "counts"))

# percentage of the total reads in airway dataset for SRR1039508 which overlaps autosome of txdb
autosome_reads/total_reads

# Output
# 0.9004193

# Question 10 
library(AnnotationHub)
ah <- AnnotationHub()
ah_E096 <- query(ah, c("E096", "H3K4me3", "narrowPeak"))
ah_record <- ah_E096[["AH30596"]]

ah_record_autosome <- keepSeqlevels(ah_record, autosome, pruning.mode = "coarse")
ah_record_ncbi <- renameSeqlevels(ah_record_autosome, txdb_ncbi)

ncbi_group <- extractSeqlevelsByGroup(species = "Homo sapiens", style = "NCBI", group = "auto")
sample_ncbi <- keepSeqlevels(range(rowRanges(sample_SRR1039508_autosome)), ncbi_group)

ov <- subsetByOverlaps(promoters(sample_ncbi), ah_record_ncbi)
ov <- subsetByOverlaps(sample_SRR1039508, ov)

median(assay(ov, "counts"))

# Output
# 205

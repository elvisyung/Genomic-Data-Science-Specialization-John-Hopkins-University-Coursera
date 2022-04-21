# Week 1 - Module 1 Quiz (Questions Relating to: GRanges, AnnotationHub and Roadmap Epigenomics)

# Question 1 
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(dbplyr)
library(BiocFileCache)
library(rtracklayer)
library(AnnotationHub)

# retrieve the cpg
ah <- AnnotationHub()
ah_human <- subset(ah, species == "Homo sapiens")
ah_human_cpg <- query(ah_human, "CpG Islands")
ah_human_cpg_record <- ah_human_cpg[["AH5086"]]

# extract autosomes
filter <- c(paste("chr", 1:22, sep=""))
split_record <- split(ah_human_cpg_record, seqnames(ah_human_cpg_record))
autosomes <- split_record[filter]

# check the number of autosomes
unlist(autosomes)

# Important Output:
# GRanges object with 26641 ranges and 1 metadata column

# Question 2
autosomes[4]

# Important Output:
# GRanges object with 1031 ranges and 1 metadata column

# Question 3
# retrieve the record
ah_H3K4me3 <- query(ah, c("H3K4me3", "narrowpeak", "E003"))

# select the narrowpeak
ah_H3K4me3_record <- ah_H3K4me3[["AH29884"]]

# extract autosomes and check the number of regions cover
split_H3K4me3 <- split(ah_H3K4me3_record, seqnames(ah_H3K4me3_record))
H3K4me3_autosomes <- split_H3K4me3[filter]
sum(width(unlist(H3K4me3_autosomes)))

# Output:
# 41135164

# Question 4
# retrieve the record
ah_H3K27me3 <- query(ah, c("H3K27me3", "narrowpeak", "E003"))
ah_H3K27me3_record <- ah_H3K27me3[["AH29892"]]

# extract autosomes
split_H3K27me3 <- split(ah_H3K27me3_record, seqnames(ah_H3K27me3_record))
H3K27me3_autosomes <- split_H3K27me3[filter]

# create a subset of extracted autosomes
ah_H3K27me3_autosomes <- subset(ah_H3K27me3_record, seqnames %in% filter)

# mean signalValue
mean_signalValue <- mean(ah_H3K27me3_autosomes$signalValue)
mean_signalValue

# Output
# 4.770728

# Question 5
# intersect between two records
bivalent <- intersect(unlist(H3K4me3_autosomes), unlist(H3K27me3_autosomes))
sum(width(bivalent))

# Output
# 10289096

# Question 6
# find bivalent regions overlap CpG Islands
cpg_autosomes <- autosomes
cpg_bivalent <- findOverlaps(bivalent, unlist(cpg_autosomes))

# calculate the fraction of the bivalent regions overlap CpG Islands
fraction_bivalent <- length(unique(queryHits(cpg_bivalent)))/length(bivalent)
fraction_bivalent

# Output
# 0.5382644

# Question 7
cpg_bivalent_intersect <- intersect(bivalent, unlist(cpg_autosomes))

# calculate the fration of the bases intersected between CpG Islands and bivalent
fraction_bivalent_intersect <- sum(width(reduce(cpg_bivalent_intersect)))/sum(width(unlist(cpg_autosomes)))
fraction_bivalent_intersect

# Output
# 0.241688

# Question 8
# extract CpG Islands within 10kb
cpg_10kb <- resize(unlist(cpg_autosomes), width = 20000 + width(unlist(cpg_autosomes)), fix = "center")

cpg_10kb_bivalent <- intersect(cpg_10kb, bivalent)
sum(width(cpg_10kb_bivalent))

# Output
# 9782086

# Question 9 
# calculate human genome size
chr_list <- c(paste("chr", 1:22, sep=""))
genome <- keepSeqlevels(ah_human_cpg_record, chr_list, pruning.mode = "coarse")
genome_size <- sum(as.numeric(seqlengths(genome)))

# calculate the fraction of human genome which contained a CpG Island
cpg_autosomes_size <- sum(as.numeric(width(unlist(cpg_autosomes))))
cpg_autosomes_size / genome_size

# Output
# 0.007047481

# Question 10
# calculate InOut matrix
overlapMat <- matrix(0,, ncol = 2, nrow = 2)
colnames(overlapMat) <- c("in", "out")
rownames(overlapMat) <- c("in", "out")
overlapMat[1,1] <- sum(width(cpg_bivalent_intersect))
overlapMat[1,2] <- sum(width(setdiff(bivalent, unlist(cpg_autosomes))))
overlapMat[2,1] <- sum(width(setdiff(unlist(cpg_autosomes), bivalent)))
overlapMat[2,2] <- genome_size - sum(overlapMat)

# calculate odds-ratio
oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2])
oddsRatio

# Output
# 169.0962






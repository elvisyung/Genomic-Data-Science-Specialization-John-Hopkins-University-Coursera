# Week 2 - Module 2 Quiz (Questions: Biostrings and BSgenomes usuage)

# Question 1
library(AnnotationHub)
library(Biostrings)
library(BSgenome)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# count total bases on chr22
alphabet_frequency <- alphabetFrequency(Hsapiens$chr22)
total_bases <- sum(alphabet_frequency[c('A','G','C','T')])

# count GC bases on chr22  
GC_bases <- sum(alphabet_frequency[c('G','C')])

#calculate GC ratio
GC_content <- GC_bases/total_bases
GC_content

# Output
# 0.4798807

# Question 2
# retrieve record
ah <- AnnotationHub()
H3K27me3_qh <- query(ah, c("H3K27me3", "E003", "narrowPeak"))
H3K27me3_record <- H3K27me3_qh[["AH29892"]]

# extract chr 22 and sequences
H3K27me3_chr22 <- subset(H3K27me3_record, seqnames == "chr22")
H3K27me3_chr22_views <- Views(Hsapiens, H3K27me3_chr22)

# calculate mean GC content
GC_contents <- letterFrequency(H3K27me3_chr22_views, "GC", as.prob = TRUE)
mean_GC <- mean(GC_contents)

mean_GC

# Output
# 0.528866

# Question 3
signal_value <- mcols(H3K27me3_chr22_views)$signalValue
cor(signal_value, GC_contents)

# Output
# 0.004467924

# Question 4
# retrieve record
H3K27me3_fc <- query(ah, c("H3K27me3", "E003", "fc.signal"))
H3K27me3_fc_record <- H3K27me3_fc[["AH32033"]]

# get subset data on chr22
gr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
H3K27me3_fc_gr <- import(H3K27me3_fc_record, which = gr22, as = "Rle")
H3K27me3_fc_gr22 <- H3K27me3_fc_gr$chr22

# view fc.signal data
fc.signal <- Views(H3K27me3_fc_gr22, start = start(H3K27me3_chr22), end = end(H3K27me3_chr22))

# calculate the correlation between the average of fc.signal and signalValue
fc.signal_mean <- mean(fc.signal)
cor(fc.signal_mean, signal_value)

# Output
# 0.9149614

# Question 5 
sum(H3K27me3_fc_gr22 >= 1)

# Ouput
# 10914671

# Question 6
H3K27me3_E055 <- query(ah, c("H3K27me3", "E055"))
H3K27me3_E055_record <- H3K27me3_E055[["AH32470"]]

# get subset data on chr22
gr_chr22 <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))
H3K27me3_fc_gr_E055 <- import(H3K27me3_E055_record, which = gr_chr22, as = "Rle")
H3K27me3_fc_gr22_E055 <- H3K27me3_fc_gr_E055$chr22

# identify region
region_E003 <- as(slice(H3K27me3_fc_gr22, upper = 0.5), "IRanges")
region_E055 <- as(slice(H3K27me3_fc_gr22_E055, lower = 2), "IRanges")
inter_region <- intersect(region_E003, region_E055)
sum(width(inter_region))

# Output
# 1869937

# Question 7
## retrieve the cpg
ah_human <- subset(ah, species == "Homo sapiens")
ah_human_cpg <- query(ah_human, "CpG Islands")
ah_human_cpg_record <- ah_human_cpg[["AH5086"]]

# get subset data on chr22
ah_human_cpg_chr22 <- subset(ah_human_cpg_record, seqnames == "chr22")
ah_human_cpg_chr22_views <- Views(Hsapiens, ah_human_cpg_chr22)

# calculate observed GC bases
observed_GC <- dinucleotideFrequency(ah_human_cpg_chr22_views)[,7]/width(ah_human_cpg_chr22_views)

# calculate expected GC bases
freq_C <- letterFrequency(ah_human_cpg_chr22_views, "C")
freq_G <- letterFrequency(ah_human_cpg_chr22_views, "G")
expected_GC <- (freq_C/width(ah_human_cpg_chr22_views))*(freq_G/width(ah_human_cpg_chr22_views))

# calculate the average observed-to-expected ratio of CpG dinucleotides
mean(observed_GC/expected_GC)

# Output
# 0.8340929

# Question 8 
TATA_boxes <- countPattern("TATAAA", Hsapiens$chr22) + countPattern("TATAAA", reverseComplement(Hsapiens$chr22))
TATA_boxes

# Output
# 27263

# Question 9
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- GRanges(seqnames = "chr22", ranges = IRanges(start = start(Hsapiens$chr22), end = end(Hsapiens$chr22)))

# find promoters of transcripts on chr 22
gr_trans_chr22 <- subsetByOverlaps(transcripts(txdb), gr, ignore.strand = TRUE)
proms <- promoters(gr_trans_chr22, upstream = 900, downstream = 100)

# find coding sequences on chr 22
gr_cds_chr22 <- subsetByOverlaps(cds(txdb), gr, ignore.strand = TRUE)

# find overlaps between promoters of transcripts and coding sequences
proms_cds <- findOverlaps(proms, gr_cds_chr22)

# calculate TATA box on overlaps
count = 0
for (i in unique(queryHits(proms_cds))){
  proms_cds_view <- Views(Hsapiens, proms[i])
  count = count + vcountPattern("TATAAA", DNAStringSet(proms_cds_view))
}

count

# Output
# 57

# Question 10
# calculate transcript lengths
trans_len_chr22 <- transcriptLengths(txdb, with.cds_len = TRUE)
trans_len_chr22 <- trans_len_chr22[trans_len_chr22$cds_len > 0,]

# find promoters from different transcripts to overlap
trans_eval <- proms[mcols(proms)$tx_id %in% trans_len_chr22$tx_id]
result = sum(coverage(trans_eval) > 1)
result["chr22"]

# Output
# 306920
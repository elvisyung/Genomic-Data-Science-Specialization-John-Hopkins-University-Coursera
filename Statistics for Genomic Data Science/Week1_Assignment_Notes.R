# Week 1 - Module 1 Quiz Notes (Questions: Exploratory Analysis and Clustering)

# Question 2
knitr::opts_chunk$set(cache=TRUE)
x = rnorm(10)
plot(x,pch=19,col="dodgerblue")
y = rbinom(20,size=1,prob=0.5)
table(y)

# Question 3
library(Biobase)
library(GenomicRanges)
data(sample.ExpressionSet, package = "Biobase")
se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)

# Question 5
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

# Question 6
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata_bm=pData(bm)

library(plotrix)
pie3D(pdata_bm$num.tech.reps,labels=pdata_bm$tissue.type)

# Question 7
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)

# Answers:
# row_sums = rowSums(edata)
# edata = edata[order(-row_sums),]
# index = 1:500
# heatmap(edata[index,],Rowv=NA,Colv=NA)

# Question 8
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)

# make an MA-plot
mm = log2(edata[,1]+1) - log2(edata[,2]+1)
aa = log2(edata[,1]+1) + log2(edata[,2]+1)
plot(aa,mm,col=2)

library(DESeq2)
rld <- rlog(exprs(bm))

y_rld = rld[,1] - rld[,2]
x_rld = rld[,1] - rld[,2]
plot(x_rld, y_rld, col = "blue", type = "p")

# Question 9
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

#With no changes to the data
dist1 = dist(t(edata))
hclust1 = hclust(dist1)

par(mar=c(0, 4, 4, 2))
plot(hclust1, hang = -1, main="origin", labels=FALSE)

#After filtering all genes with rowMeans less than 100
low_genes = rowMeans(edata) < 100
filter_edata = filter(edata, !low_genes)
f_dist1 = dist(t(filter_edata))
f_hclust1 = hclust(f_dist1)

par(mar=c(0, 4, 4, 2))
plot(f_hclust1, hang = -1, main="remove low expression", labels=FALSE)

#After taking the log2 transform of the data without filtering
log_edata = log2(edata + 1)
l_dist1 = dist(t(log_edata))
l_hclust1 = hclust(l_dist1)

par(mar=c(0, 4, 4, 2))
plot(l_hclust1, hang=-1, main="perform log2 transform", labels=FALSE)

# Question 10 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata = log2(edata + 1)

# perfrom k-means clustering
set.seed(1235)
k2 = kmeans(edata,centers=2)
matplot(t(k2$centers),col=1:2,type="l",lwd=3)

dist1 = dist(t(edata))
hclust1 = hclust(dist1)
tree = cutree(hclust1, 2)

par(mar=c(0, 4, 4, 2))
plot(hclust1, tree, main="cutree")





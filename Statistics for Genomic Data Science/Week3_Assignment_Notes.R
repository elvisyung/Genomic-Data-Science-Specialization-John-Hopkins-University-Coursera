# Week 3 - Module 3 Quiz Notes (Questions: Logistic regressiona and multiple testing)

# Question 1
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# Fit a linear model and a logistic regression model to the data for the 3rd SNP. What are the coefficients for the SNP variable?
# recode 0 values to NA
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA

# fit a linear model
lm3 = lm(status ~ snp3)
tidy(lm3)

# Output:
# # A tibble: 2 x 5
#   term        estimate std.error statistic  p.value
#   <chr>          <dbl>     <dbl>     <dbl>    <dbl>
# 1 (Intercept)   0.544     0.0549     9.91  3.75e-22
# 2 snp3         -0.0394    0.0468    -0.842 4.00e- 1

# fit a logistic regression model
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)

# Output:
# # A tibble: 2 x 5
#   term        estimate std.error statistic p.value
#   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
# 1 (Intercept)    0.177     0.220     0.806   0.420
# 2 snp3          -0.158     0.188    -0.841   0.400

# Question 2
par(mfrow=c(1,2))

plot(status ~ snp3,pch=19)
abline(lm3,col="darkgrey",lwd=5)
plot(glm3$residuals)

# Question 3
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# fit a logistic regression model
snp10 = as.numeric(snpdata[,10])
snp10[snp10==0] = NA
glm10 = glm(status ~ snp10, family="binomial")
tidy(glm10)

# Output:
# # A tibble: 2 x 5
#   term        estimate std.error statistic p.value
#   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
# 1 (Intercept) -0.00751    0.174    -0.0433   0.965
# 2 snp10        0.00201    0.0933    0.0215   0.983

snp10_dom = (snp10 == 2)
glm10_dom = glm(status ~ snp10_dom, family="binomial")
tidy(glm10_dom)

# Output:
# # A tibble: 2 x 5
#   term          estimate std.error statistic p.value
#   <chr>            <dbl>     <dbl>     <dbl>   <dbl>
# 1 (Intercept)    -0.0339    0.0868    -0.391   0.696
# 2 snp10_domTRUE   0.0643    0.127      0.505   0.614

# Question 4
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# fit an additive logistic regression model to each SNP
results = rep(NA, dim(snpdata)[2])
for (i in 1:ncol(snpdata)){
  snpdata_i = as.numeric(snpdata[,i])
  snpdata_i[snpdata_i == 0] = NA
  glm_i = glm(status ~ snpdata_i, family = "binomial")
  results[i] = tidy(glm_i)$statistic[2]
}

# average effect size
mean(results)

# 0.007155377

# minimum effect size
min(results)

# -4.251469

# maximum effect size
max(results)

# 3.900891

# Question 5
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# square the coefficients
results_coeff_squre =  results^2

# correlation with the results from using snp.rhs.tests and chi.squared
glm_all = snp.rhs.tests(status ~ 1, snp.data = sub.10)
cor(results_coeff_squre, chi.squared(glm_all))

# 0.9992946

# Question 6
library(ballgown)
library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# Do the log2(data + 1) transform and fit calculate F-statistics for the difference between studies/populations using genefilter:rowFtests and using genefilter:rowttests.
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

edata = log2(as.matrix(edata) + 1)

# perform rowttests
tstats_obj = rowttests(edata, as.factor(pdata$population))
tidy(tstats_obj)

# Output:
# # A tibble: 3 x 13
#   column     n   mean    sd   median trimmed     mad       min   max range
#   <chr>  <dbl>  <dbl> <dbl>    <dbl>   <dbl>   <dbl>     <dbl> <dbl> <dbl>
# 1 stati… 12984 -3.43  3.88  -2.84    -3.16   2.94    -2.31e+ 1  8.31 31.5 
# 2 dm     52580 -0.129 0.372  0       -0.0283 0       -4.37e+ 0  1.51  5.89
# 3 p.val… 12984  0.168 0.266  0.00441  0.109  0.00441  2.23e-47  1.00  1.00
# # … with 3 more variables: skew <dbl>, kurtosis <dbl>, se <dbl>

# perform rowFtests
fstats_obj = rowFtests(edata, as.factor(pdata$population))
tidy(fstats_obj)

# # A tibble: 2 x 13
#   column     n   mean     sd  median trimmed     mad      min    max  range
#   <chr>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>    <dbl>  <dbl>  <dbl>
# 1 stati… 12984 26.8   39.5   8.40     18.4   8.20    3.39e- 7 535.   535.  
# 2 p.val… 12984  0.168  0.266 0.00441   0.109 0.00441 2.23e-47   1.00   1.00
# # … with 3 more variables: skew <dbl>, kurtosis <dbl>, se <dbl>

par(mfrow=c(1,2))
hist(tstats_obj$statistic, col=2)
hist(fstats_obj$statistic, col=2)

# Question 7
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

# First test for differences between the studies using the DESeq2 package using the DESeq function. Then do the log2(data + 1) transform and do the test for differences between studies using the limma package and the lmFit, ebayes and topTable functions.
library(DESeq2)
library(limma)
library(edge)
library(genefilter)

# using DESeq2 test the differences between the studies
de = DESeqDataSetFromMatrix(edata, pdata, ~study)
glm_de = DESeq(de)
result_de = results(glm_de)

# using limma test the differences
edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ as.factor(pdata$study))
fit_limma = lmFit(edata, mod)
ebayes_limma = eBayes(fit_limma) 
top = topTable(ebayes_limma,number=dim(edata)[1], sort.by="none")

# correlation in the statistics between two analyses
cor(result_de$stat, top$t)

# 0.9278568

# make an MA-plot
y = cbind(result_de$stat, top$t)
limma::plotMA(y)

# Question 8
# DESeq analysis
fp_bh = p.adjust(result_de$pvalue, method="BH")
sum(fp_bh < 0.05)

# 1995

# limma analysis
fp_bh = p.adjust(top$P.Value, method="BH")
sum(fp_bh < 0.05)

# 2807




# Week 2 - Module 2 Quiz Notes (Questions: Pre-processing, linear regression and batch effects)

# Question 1
library(ballgown)
library(Biobase)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# No transformations
svd1 = svd(edata)
ori_pca = svd1$d^2/sum(svd1$d^2)
ori_pca[1]

# 0.8873421

# log2 transform
edata_log2 = log2(edata + 1)
svd2 = svd(edata_log2)
log2_pca = svd2$d^2/sum(svd2$d^2)
log2_pca[1]

# 0.9737781

# log2 transform, subtract row means
edata_centered = edata_log2 - rowMeans(edata_log2)
svd3 = svd(edata_centered)
centered_data_pca = svd3$d^2/sum(svd3$d^2)
centered_data_pca[1]

# 0.3463729

# Question 2
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# Performing the log2(data + 1)transform and subtract row means from the samples. 
edata_log2 = log2(edata + 1)
edata_centered = edata_log2 - rowMeans(edata_log2)

# use svd to calculate the singular vectors
set.seed(333)
svd1 = svd(edata_centered)

edata_kmeans = kmeans(t(edata_centered), centers=2)
cor.test(svd1$v[,1], edata_kmeans$cluster)

# Output:
# Pearson's product-moment correlation
# 
# data:  svd1$v[, 1] and edata_kmeans$cluster
# t = -19.683, df = 127, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9049326 -0.8176191
# sample estimates:
#        cor 
# -0.8678247

# Question 3
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# fit linear model
lm1 = lm(edata[1,] ~ pdata_bm$num.tech.reps)

# plot the data
plot(pdata_bm$num.tech.reps,edata[1,])
abline(lm1$coeff[1], lm1$coeff[2], col=2, lwd=3)

# Question 4
# fit linear model
lm2 = lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
summary(lm2)

# Output:
# Call:
# lm(formula = edata[1, ] ~ pdata_bm$age + pdata_bm$gender)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -734.35 -229.31   -3.26  243.02  768.09 
# 
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      2331.581    438.181   5.321 0.000139 ***
# pdata_bm$age      -23.913      6.488  -3.686 0.002744 ** 
# pdata_bm$genderM -207.257    236.431  -0.877 0.396610    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 469.6 on 13 degrees of freedom
#   (3 observations deleted due to missingness)
# Multiple R-squared:  0.5147, Adjusted R-squared:   0.44 
# F-statistic: 6.894 on 2 and 13 DF,  p-value: 0.009102

# Question 5
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

# Perform the log2(data + 1) transformation 
edata = log2(edata + 1)

# fit a regression model to each sample, using population as the outcome
mod = model.matrix(~ pdata$population)
fit = lm.fit(mod, t(edata))

# dimension of the residual matrix
dim(fit$residuals)

# 129 52580

# dimension of the effects matrix
dim(fit$effects)

# 129 52580

# dimension of the coefficients matrix
dim(fit$coefficients)

# 2 52580

# Question 7
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# Fit regression models to the expresison data where age is the outcome variable using the Lmfit function from the limma package.
library(devtools)
library(Biobase)
library(limma)
library(edge)

# subset the expression data to the samples without mimssing values of age
pdata_bm = na.omit(pdata_bm)
edata = edata[,rownames(pdata_bm), drop=FALSE]

# fit many regression models to the expression data where age is the outcome
mod_adj = model.matrix(~ pdata_bm$age)
fit_limma = lmFit(edata,mod_adj)

fit_limma$coefficients[1000,]

# Output:
#  (Intercept) pdata_bm$age 
#   2469.87375    -27.61178

# make a plot of the 1,000th gene and fitted values
intercept = fit_limma$coefficients[1000,][1]
slope = fit_limma$coefficients[1000,][2]
x = edata[1000,]*slope+intercept

plot(x,pdata_bm$age)

# Question 8 
pdata_bm$tissue.type

# Output:
#  [1] adipose          adrenal          brain            breast          
#  [5] colon            heart            kidney           liver           
#  [9] lung             lymphnode        ovary            prostate        
# [13] skeletal_muscle  testes           thyroid          white_blood_cell
# 17 Levels: adipose adrenal brain breast colon heart kidney liver ... white_blood_cell

# Question 9
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

head(pdata)

# Output:
#         sample.id num.tech.reps population      study
# NA06985   NA06985             1        CEU Montgomery
# NA06986   NA06986             1        CEU Montgomery
# NA06994   NA06994             1        CEU Montgomery
# NA07000   NA07000             1        CEU Montgomery
# NA07037   NA07037             1        CEU Montgomery
# NA07051   NA07051             1        CEU Montgomery

# Question 10
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

# Using command set.sedd(33353) to estimate a single surrogate variable using the sva function after log2(data + 1) transforming the expression data, removing rows with rowMeans less than 1, and treating age as the outcome.
library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

# preprocessing the data
set.seed(33353)
pheno = na.omit(pdata_bm)
edata = edata[,rownames(pheno), drop=FALSE]
edata = log2(edata + 1)
edata = edata[rowMeans(edata) > 1,]

# fit a sva model
mod = model.matrix(~age, data=pheno)
mod0 = model.matrix(~1, data=pheno)
sva1 = sva(edata, mod,mod0, n.sv=2)

# Output:
# Number of significant surrogate variables is:  2 
# Iteration (out of 5 ):1  2  3  4  5

# correlation between surrogate for batch and age
cor(sva1$sv, pheno$age)

#            [,1]
# [1,] -0.1965417
# [2,] -0.1560322

# correlation between surrogate for batch and race
cor(sva1$sv, as.numeric(pheno$race))

#             [,1]
# [1,] -0.35780610
# [2,]  0.04154497





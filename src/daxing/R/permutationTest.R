library("ggplot2"); theme_set(theme_bw())
library("lmPerm")
library("coin")
library("gtools")
setwd("/Users/xudaxing/Documents/deleteriousMutation/001_vmap2.1Before20200525/003_rScript/001_data")
ants <- read.table("permutationTest.txt",header = T)
print(ants)
#plot
col<-c("red","blue")
p <- ggplot(ants,aes(x=place, y= colonies))+
  geom_boxplot(position = position_dodge(0.8),outlier.colour = 'black',alpha=0.8)+
  stat_boxplot(geom = "errorbar",width=0.12,position = position_dodge(0.5))+
  stat_sum()+
  geom_point(position = position_jitterdodge(0.8),alpha = 0.6,aes(col=place),size=0.3)+
  theme(panel.background = element_rect(size = 0.7, colour = 'black',fill = 'white'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),text = element_text(size = 12),legend.position = 'none')+
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(size=0.4, colour = "black"),text = element_text(size = 12))+
  # xlab("")+ylab("Relative contribution (%)")+
  scale_fill_manual(values = col)+
  scale_color_manual(values = col)+
  # ggtitle("M111")+ 
  theme(plot.title = element_text(hjust = 0.5)) 
p

##########################################################################################################################
########################   permutation significance test by mean    ####################################################
##########################################################################################################################
setwd("/Users/xudaxing/Documents/deleteriousMutation/001_vmap2.1Before20200525/003_rScript/001_data")
ants <- read.table("permutationTest.txt",header = T)
set.seed(101) ## for reproducibility
nsim <- 1000
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(ants)) ## by default, 不放回抽样
  bdat <- transform(ants,colonies=colonies[perm])
  ## compute & store difference in means; store the value
  res[i] <- mean(bdat[bdat$place=="field","colonies"])-
    mean(bdat[bdat$place=="forest","colonies"])
}
obs <- mean(ants[ants$place=="field","colonies"])-
  mean(ants[ants$place=="forest","colonies"])
## append the observed value to the list of results
# res <- c(res,obs)
# hist(res,col="gray",las=1,main="")
# abline(v=obs,col="red")
#permutation res plot
par(las=1,bty="l")
plot(prop.table(table(round(res,2))),
     ylab="Proportion",axes=FALSE)
axis(side=2)
points(obs,0,pch=16,cex=1.5,col="red")
######   permutation significance test by mean
2*mean(res>=obs)          ## doubling (as suggested by JD)        double tails
mean(abs(res)>=abs(obs))  ## count both tails: matches lmPerm     double tails

##########################################################################################################################
########################   permutation significance test by t test    ####################################################
##########################################################################################################################
ants <- read.table("permutationTest.txt",header = T)
tt <- t.test(colonies~place,data=ants,var.equal=TRUE)
obs<-tt$statistic
set.seed(101) ## for reproducibility
nsim <- 1000
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(ants)) ## by default, 不放回抽样
  bdat <- transform(ants,colonies=colonies[perm])
  tt <- t.test(colonies~place, data=bdat, var.equal = T)
  res[i] <- tt$statistic
}
#plot null distribution
par(las=1,bty="l")
plot(prop.table(table(round(res,2))),
     ylab="Proportion",axes=FALSE)
axis(side=1)
points(obs,0,pch=16,cex=1.5,col="red")
# permutation significance
res <- c(res,obs)
2*mean(res>=obs)          ## doubling (as suggested by JD)       double tails
mean(abs(res)>=abs(obs))  ## count both tails: matches lmPerm    double tails


##########################################################################################################################
########################   permutation significance test by p value    ####################################################
##########################################################################################################################
ants <- read.table("permutationTest.txt",header = T)
tt <- t.test(colonies~place,data=ants,var.equal=TRUE)
obs<-tt$p.value
set.seed(101) ## for reproducibility
nsim <- 10000
res <- numeric(nsim) ## set aside space for results
for (i in 1:nsim) {
  ## standard approach: scramble response value
  perm <- sample(nrow(ants)) ## by default, 不放回抽样
  bdat <- transform(ants,colonies=colonies[perm])
  tt <- t.test(colonies~place, data=bdat, var.equal = T)
  res[i] <- tt$p.value
}
res <- c(res,obs)
#plot null distribution
par(las=1,bty="l")
plot(prop.table(table(round(res,2))),
     ylab="Proportion",axes=FALSE)
axis(side=2)
points(obs,0,pch=16,cex=1.5,col="red")
# permutation significance
2*mean(res<=obs)          ## doubling (as suggested by JD)        double tails
mean(abs(res)<=abs(obs))  ## count both tails: matches lmPerm     double tails


##########################################################################################################################
########################   permutation significance test by coin package    ####################################################
##########################################################################################################################
library(survival)
library(coin)
setwd("/Users/xudaxing/Documents/deleteriousMutation/001_vmap2.1Before20200525/003_rScript/001_data")
table <- read.table("permutationTest.txt",header = T,stringsAsFactors = T)
## default: asymptotic
table$place = factor(table$place,levels = c("forest","field"))
levels(table$place)
oneway_test(colonies~place,data=table,alternative = "less")
## exact distribution
oneway_test(colonies~place,data=table,distribution="exact")
## approximate (random-sampling/Monte Carlo)
oneway_test(colonies~place,data=table,distribution=approximate(10000))

##########################################################################################################################
########################   Kruskal-Wallis test by rank     ####################################################
##########################################################################################################################
library(dplyr)
library(ggplot2)
setwd("/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/003_permutationTest")
table<-read.table("deleteriousLoad.txt",header = T,stringsAsFactors = T)
head(table)
levels(table$Subspecies)
table$Subspecies <- ordered(table$Subspecies, levels = c("Wild_emmer", "Domesticated_emmer", "Free_threshing_tetraploid","Landrace","Cultivar","Ae.tauschii"))
col<-c("#ffd702","#7f5701","#016699","#fc6e6e","#9900ff","#87cef9")

#load 回归校正前
p<-ggplot(table,aes(x=Sub,y=Ratio_delVSsyn,fill=Subspecies))
p+geom_boxplot()+scale_fill_manual(values = col)+
  labs(x="Sub",y="Del/syn per accession")+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5),legend.position = "NA")
# kruskal.test one factor
kruskal.test(Ratio_delVSsyn ~ Subspecies, data = table)
pairwise.wilcox.test(table$Ratio_delVSsyn, table$Subspecies,p.adjust.method = "BH")
# two-way ANOVA test
res.aov2 <- aov(Ratio_delVSsyn ~ Sub * Subspecies, data = table)
summary(res.aov2)
#permutation test
levels(table$Sub)
oneway_test(Ratio_delVSsyn~Sub,data=table)

#load 回归校正后
p<-ggplot(table,aes(x=Sub,y=resi,fill=Subspecies))
p+geom_boxplot()+scale_fill_manual(values = col)+
  labs(x="Sub",y="Corrected del/syn per accession")+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5),legend.position = "NA")
# kruskal.test one factor
kruskal.test(resi ~ Subspecies, data = table)
pairwise.wilcox.test(table$resi, table$Subspecies,p.adjust.method = "BH")
# two-way ANOVA test
res.aov2 <- aov(resi ~ Sub * Subspecies, data = table)
summary(res.aov2)
#permutation test
oneway_test(Ratio_delVSsyn~Sub,data=table,alternative = "less")



##########################################################################################################################
########################   pairwisePermutationTest and pairwisePermutationMatrix   ####################################################
##########################################################################################################################

### One-way permutation test, hypothetical data
Input =("
Factor Response
  A     4.6
  A     5.5
  A     3.4
  A     5.0
  A     3.9
  A     4.5
  B     3.6
  B     4.5
  B     2.4
  B     4.0
  B     2.9
  B     3.5
  C     2.6
  C     3.5
  C     1.4
  C     3.0
  C     1.9
  C     2.5
  D     4.7
  D     5.6
  D     3.5
  D     5.1
  D     4.0
  D     4.6
")

Data = read.table(textConnection(Input),header=TRUE)


### Order groups by median; pairwisePermutationTest
Data$Factor = factor(Data$Factor,levels = c("D", "A", "B", "C"))
headtail(tableA)
PT = pairwisePermutationTest(Response ~ Factor, data = Data, method="fdr")
PT
cldList(p.adjust ~ Comparison,data = PT,threshold  = 0.05)
### Order groups by median; pairwisePermutationMatrix
Data$Factor = factor(Data$Factor,levels = c("D", "A", "B", "C"))
PM = pairwisePermutationMatrix(Response ~ Factor, data = Data, method="fdr")
PM
### multcompLetters
multcompLetters(PM$Adjusted,compare="<",threshold=0.05,Letters=letters,reversed = FALSE)

#####################################################################################################################################3
library(tidyverse)
library(caret)
theme_set(theme_bw())
# Load the data
data("marketing", package = "datarium")
# Inspect the data
sample_n(marketing, 3)
# Split the data into training and test set
set.seed(123)
training.samples <- marketing$sales %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data  <- marketing[training.samples, ]
test.data <- marketing[-training.samples, ]


# Build the model
model <- lm(sales ~., data = train.data)
# Summarize the model
summary(model)
# Make predictions
predictions <- model %>% predict(test.data)
# Model performance
# (a) Prediction error, RMSE
RMSE(predictions, test.data$sales)
# (b) R-square
R2(predictions, test.data$sales)













library(survival)
library(coin)
library(FSA)
library(rcompanion)
library(multcompView)
library(ggplot2)
library(dplyr)
### --------------------------------------------------------------
setwd("/Users/xudaxing/Documents/deleteriousMutation/003_vmap2.1_20200628/004_deleterious/003_permutationTest")
table <- read.table("Correct_resi_byIBS_catogery.txt",header = T,stringsAsFactors = F)
table$Subspecies <- ordered(table$Subspecies, levels = c("Ae.tauschii","Wild_emmer", "Domesticated_emmer", "Free_threshing_tetraploid","Landrace","Cultivar"))
col<-c("#87cef9","#ffd702","#7f5701","#016699","#fc6e6e","#9900ff")

#回归校正前 del load
p<-ggplot(table,aes(x=Sub,y=Ratio_delVSsyn,fill=Subspecies))
p+geom_boxplot()+scale_fill_manual(values = col)+
  labs(x="Sub",y="Del/syn per accession")+
  theme(panel.background = element_blank(),axis.line = element_line(size = 0.5),legend.position = "NA",text = element_text(size = 20))
#回归校正后 del load
p<-ggplot(table,aes(x=Sub,y=resi,fill=Subspecies))
q<-p+geom_boxplot(position = position_dodge(0.8))+scale_fill_manual(values = col)+
  labs(x="Sub",y="Del/syn per accession")+
  scale_y_continuous(limits = c(-0.0025,0.027))+
  # geom_text(data = df,aes(x=Sub,fill))
  theme_bw()+
  theme(panel.background = element_blank(),panel.border = element_blank(),
        axis.line = element_line(size=0.5),legend.position = "NA",
        text = element_text(size = 20),panel.grid = element_blank())
q
#

### new df
### summary
library(dplyr)
df_group <- group_by(table,Sub,Subspecies)
out <- dplyr::summarise(df_group,mean = mean(Ratio_delVSsyn), sd = sd(Ratio_delVSsyn))
signifi <- c("A","B","B","BC","C","A","B","B","B","B","A","B","B")
permu <- c("a","a","b","b","c","a","b","b","b","b","a","b","b")
out <- mutate(out,y1=-0.001,y2=-0.0025)
out$signifi <- signifi
out$permu <- permu
out
####### final plot on correction boxplot with significant character
q <- p + geom_text(data=out,aes(x=Sub,y=y1,group=Subspecies,label=signifi),position = position_dodge(0.8),size=2.5)+
  geom_text(data=out,aes(x=Sub,y=y2,group=Subspecies,label=permu),position = position_dodge(0.8),size=3)
q
ggsave("correctedDelLoad.pdf",width = 5,height = 5)

##############======================================================================================================
# 数据清洗 tidyr
# 长表 宽表 
# summary: dpylr
# df <- group_by(table,Sub,Subspecies)
out <- dplyr::summarise(df_group,r = cor(Ratio_delVSsyn,IBSdistancetoCS2017))




out
#############################################################################################
#-----------------------------------    回归校正前 del load  ----------------------------
#############################################################################################

###------------ A sub
tableA<-filter(table, Sub=="A")
tableA$Subspecies <- factor(tableA$Subspecies, levels = c("Wild_emmer", "Domesticated_emmer", "Free_threshing_tetraploid","Landrace","Cultivar"))
#pairwise.wilcox.test  
wtAtwo.sided<-pairwise.wilcox.test(tableA$Ratio_delVSsyn, tableA$Subspecies,p.adjust.method = "BH",alternative="two.sided")
wtAtwo.sided$p.value/2
### One-way permutation test
independence_test(tableA$Ratio_delVSsyn ~ tableA$Subspecies,data = tableA)
pmAtwo<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Subspecies, data=tableA, method = "fdr", alternative="two.sided")
pmAtwo$Adjusted

library(corrplot)
library(RColorBrewer)
#order="hclust"
corrplot(pmAtwo$Adjusted, method = "number",type="upper", is.corr = F, col=brewer.pal(n=8, name="PuOr"))


###------------ B sub
tableB<-filter(table, Sub=="B")
#多重比较前要按按中位数排序
tableB$Subspecies <- factor(tableB$Subspecies, levels = c("Wild_emmer", "Domesticated_emmer", "Free_threshing_tetraploid","Landrace","Cultivar"))
#pairwise.wilcox.test
wtBtwo<-pairwise.wilcox.test(tableB$Ratio_delVSsyn, tableB$Subspecies,p.adjust.method = "BH",alternative="two.sided")
wtBtwo
### One-way permutation test
independence_test(tableB$Ratio_delVSsyn ~ tableB$Subspecies,data = tableB)
pmBtwo<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Subspecies, data=tableB, method = "fdr",alternative="two.sided")
pmBtwo

###------------ D sub
tableD<-filter(table, Sub=="D")
#多重比较前要按按中位数排序
tableD$Subspecies <- factor(tableD$Subspecies, levels = c("Ae.tauschii","Landrace","Cultivar"))
#pairwise.wilcox.test
wtDtwo<-pairwise.wilcox.test(tableD$Ratio_delVSsyn, tableD$Subspecies,p.adjust.method = "BH",alternative="two.sided")
wtDtwo
### One-way permutation test
independence_test(tableD$Ratio_delVSsyn ~ tableD$Subspecies,data = tableD)
pmDtwo<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Subspecies, data=tableD, method = "fdr",alternative="two.sided")
pmDtwo

###------------ Wild_emmer A B
tableWildEmmer<-filter(table,table$Subspecies=="Wild_emmer")
#多重比较前要按按中位数排序
tableWildEmmer$Sub <- factor(tableWildEmmer$Sub, levels = c("A","B"))
#pairwise.wilcox.test
wtWildEmmer<-pairwise.wilcox.test(tableWildEmmer$Ratio_delVSsyn, tableWildEmmer$Sub,p.adjust.method = "BH",alternative="two.sided")
wtWildEmmer
### One-way permutation test
independence_test(tableWildEmmer$Ratio_delVSsyn ~ tableWildEmmer$Sub,data = tableWildEmmer)
pmWildEmmer<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Sub, data=tableWildEmmer, method = "fdr",alternative="two.sided")
pmWildEmmer

###------------ Domesticated_emmer A B
tableDomesticated_emmer<-filter(table,table$Subspecies=="Domesticated_emmer")
#多重比较前要按按中位数排序
tableDomesticated_emmer$Sub <- factor(tableDomesticated_emmer$Sub, levels = c("A","B"))
#pairwise.wilcox.test
wtDomesticated_emmer<-pairwise.wilcox.test(tableDomesticated_emmer$Ratio_delVSsyn, tableDomesticated_emmer$Sub,p.adjust.method = "BH",alternative="two.sided")
wtDomesticated_emmer
### One-way permutation test
independence_test(tableDomesticated_emmer$Ratio_delVSsyn ~ tableDomesticated_emmer$Sub,data = tableDomesticated_emmer)
pmDomesticated_emmer<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Sub, data=tableDomesticated_emmer, method = "fdr",alternative="two.sided")
pmDomesticated_emmer

###------------ Free_threshing_tetraploid A B
tableFree_threshing_tetraploid<-filter(table,table$Subspecies=="Free_threshing_tetraploid")
#多重比较前要按按中位数排序
tableFree_threshing_tetraploid$Sub <- factor(tableFree_threshing_tetraploid$Sub, levels = c("A","B"))
#pairwise.wilcox.test
wtFree_threshing_tetraploid<-pairwise.wilcox.test(tableFree_threshing_tetraploid$Ratio_delVSsyn, tableFree_threshing_tetraploid$Sub,p.adjust.method = "BH",alternative="two.sided")
wtFree_threshing_tetraploid
### One-way permutation test
independence_test(tableFree_threshing_tetraploid$Ratio_delVSsyn ~ tableFree_threshing_tetraploid$Sub,data = tableFree_threshing_tetraploid)
pmFree_threshing_tetraploid<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Sub, data=tableFree_threshing_tetraploid, method = "fdr",alternative="two.sided")
pmFree_threshing_tetraploid

###------------ Landrace A B D
tableLandrace<-filter(table,table$Subspecies=="Landrace")
#多重比较前要按按中位数排序
tableLandrace$Sub <- factor(tableLandrace$Sub, levels = c("A","B","D"))
#pairwise.wilcox.test
wtLandrace<-pairwise.wilcox.test(tableLandrace$Ratio_delVSsyn, tableLandrace$Sub,p.adjust.method = "BH",alternative="two.sided")
wtLandrace
### One-way permutation test
independence_test(tableLandrace$Ratio_delVSsyn ~ tableLandrace$Sub,data = tableLandrace)
pmLandrace<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Sub, data=tableLandrace, method = "fdr",alternative="two.sided")
pmLandrace

###------------ Cultivar A B D
tableCultivar<-filter(table,table$Subspecies=="Cultivar")
#多重比较前要按按中位数排序
tableCultivar$Sub <- factor(tableCultivar$Sub, levels = c("A","B","D"))
#pairwise.wilcox.test
wtCultivar<-pairwise.wilcox.test(tableCultivar$Ratio_delVSsyn, tableCultivar$Sub,p.adjust.method = "BH",alternative="two.sided")
wtCultivar
### One-way permutation test
independence_test(tableCultivar$Ratio_delVSsyn ~ tableCultivar$Sub,data = tableCultivar)
pmCultivar<-pairwisePermutationMatrix(Ratio_delVSsyn ~ Sub, data=tableCultivar, method = "fdr",alternative="two.sided")
pmCultivar


#############################################################################################
#-----------------------------------    回归校正后 del load  ----------------------------
#############################################################################################

###------------ A sub
tableA<-filter(table, Sub=="A")
#多重比较前要按按中位数排序
tableA$Subspecies <- factor(tableA$Subspecies, levels = c("Wild_emmer", "Domesticated_emmer", "Free_threshing_tetraploid","Landrace","Cultivar"))
#pairwise.wilcox.test  
wtAtwo.sided<-pairwise.wilcox.test(tableA$resi, tableA$Subspecies,p.adjust.method = "BH",alternative="two.sided")
wtAtwo.sided$p.value/2

### One-way permutation test
independence_test(tableA$resi ~ tableA$Subspecies,data = tableA)
pmAtwo<-pairwisePermutationMatrix(resi ~ Subspecies, data=tableA, method = "fdr", alternative="two.sided")
pmAtwo
pmAtwo$Adjusted/2
multcompLetters(pmAtwo$Adjusted/2,compare="<",threshold=0.05,Letters=letters,reversed = FALSE)
###------------ B sub
tableB<-filter(table, Sub=="B")
#多重比较前要按按中位数排序
tableB$Subspecies <- factor(tableB$Subspecies, levels = c("Wild_emmer", "Domesticated_emmer", "Free_threshing_tetraploid","Landrace","Cultivar"))
#pairwise.wilcox.test
wtBtwo<-pairwise.wilcox.test(tableB$resi, tableB$Subspecies,p.adjust.method = "BH",alternative="two.sided")
wtBtwo$p.value/2

### One-way permutation test
independence_test(tableB$resi ~ tableB$Subspecies,data = tableB)
pmBtwo<-pairwisePermutationMatrix(resi ~ Subspecies, data=tableB, method = "fdr",alternative="two.sided")
pmBtwo$Adjusted/2
multcompLetters(pmBtwo$Adjusted/2,compare="<",threshold=0.05,Letters=letters,reversed = FALSE)
###------------ D sub
tableD<-filter(table, Sub=="D")
#多重比较前要按按中位数排序
tableD$Subspecies <- factor(tableD$Subspecies, levels = c("Ae.tauschii","Landrace","Cultivar"))
#pairwise.wilcox.test
wtDtwo<-pairwise.wilcox.test(tableD$resi, tableD$Subspecies,p.adjust.method = "BH",alternative="two.sided")
wtDtwo$p.value/2
### One-way permutation test
independence_test(tableD$resi ~ tableD$Subspecies,data = tableD)
pmDtwo<-pairwisePermutationMatrix(resi ~ Subspecies, data=tableD, method = "fdr",alternative="two.sided")
pmDtwo$Adjusted/2
multcompLetters(pmDtwo$Adjusted/2,compare="<",threshold=0.05,Letters=letters,reversed = FALSE)


###------------ Wild_emmer A B
tableWildEmmer<-filter(table,table$Subspecies=="Wild_emmer")
#多重比较前要按按中位数排序
tableWildEmmer$Sub <- factor(tableWildEmmer$Sub, levels = c("A","B"))
#pairwise.wilcox.test
wtWildEmmer<-pairwise.wilcox.test(tableWildEmmer$resi, tableWildEmmer$Sub,p.adjust.method = "BH",alternative="two.sided")
wtWildEmmer
### One-way permutation test
independence_test(tableWildEmmer$resi ~ tableWildEmmer$Sub,data = tableWildEmmer)
pmWildEmmer<-pairwisePermutationMatrix(resi ~ Sub, data=tableWildEmmer, method = "fdr",alternative="two.sided")
pmWildEmmer$Adjusted/2

###------------ Domesticated_emmer A B
tableDomesticated_emmer<-filter(table,table$Subspecies=="Domesticated_emmer")
#多重比较前要按按中位数排序
tableDomesticated_emmer$Sub <- factor(tableDomesticated_emmer$Sub, levels = c("A","B"))
#pairwise.wilcox.test
wtDomesticated_emmer<-pairwise.wilcox.test(tableDomesticated_emmer$resi, tableDomesticated_emmer$Sub,p.adjust.method = "BH",alternative="two.sided")
wtDomesticated_emmer
### One-way permutation test
independence_test(tableDomesticated_emmer$resi ~ tableDomesticated_emmer$Sub,data = tableDomesticated_emmer)
pmDomesticated_emmer<-pairwisePermutationMatrix(resi ~ Sub, data=tableDomesticated_emmer, method = "fdr",alternative="two.sided")
pmDomesticated_emmer$Adjusted/2

###------------ Free_threshing_tetraploid A B
tableFree_threshing_tetraploid<-filter(table,table$Subspecies=="Free_threshing_tetraploid")
#多重比较前要按按中位数排序
tableFree_threshing_tetraploid$Sub <- factor(tableFree_threshing_tetraploid$Sub, levels = c("A","B"))
#pairwise.wilcox.test
wtFree_threshing_tetraploid<-pairwise.wilcox.test(tableFree_threshing_tetraploid$resi, tableFree_threshing_tetraploid$Sub,p.adjust.method = "BH",alternative="two.sided")
wtFree_threshing_tetraploid
### One-way permutation test
independence_test(tableFree_threshing_tetraploid$resi ~ tableFree_threshing_tetraploid$Sub,data = tableFree_threshing_tetraploid)
pmFree_threshing_tetraploid<-pairwisePermutationMatrix(resi ~ Sub, data=tableFree_threshing_tetraploid, method = "fdr",alternative="two.sided")
pmFree_threshing_tetraploid$Adjusted/2
  

###------------ Landrace A B D
tableLandrace<-filter(table,table$Subspecies=="Landrace")
#多重比较前要按按中位数排序
tableLandrace$Sub <- factor(tableLandrace$Sub, levels = c("A","B","D"))
#pairwise.wilcox.test
wtLandrace<-pairwise.wilcox.test(tableLandrace$resi, tableLandrace$Sub,p.adjust.method = "BH",alternative="two.sided")
wtLandrace
### One-way permutation test
independence_test(tableLandrace$resi ~ tableLandrace$Sub,data = tableLandrace)
pmLandrace<-pairwisePermutationMatrix(resi ~ Sub, data=tableLandrace, method = "fdr",alternative="two.sided")
pmLandrace$Adjusted/2

###------------ Cultivar A B D
tableCultivar<-filter(table,table$Subspecies=="Cultivar")
#多重比较前要按按中位数排序
tableCultivar$Sub <- factor(tableCultivar$Sub, levels = c("A","B","D"))
#pairwise.wilcox.test
wtCultivar<-pairwise.wilcox.test(tableCultivar$resi, tableCultivar$Sub,p.adjust.method = "BH",alternative="two.sided")
wtCultivar
### One-way permutation test
independence_test(tableCultivar$resi ~ tableCultivar$Sub,data = tableCultivar)
pmCultivar<-pairwisePermutationMatrix(resi ~ Sub, data=tableCultivar, method = "fdr",alternative="two.sided")
pmCultivar$Adjusted/2





# x is a matrix containing the data
# method : correlation method. "pearson"" or "spearman"" is supported
# removeTriangle : remove upper or lower triangle
# results :  if "html" or "latex"
# the results will be displayed in html or latex format
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 






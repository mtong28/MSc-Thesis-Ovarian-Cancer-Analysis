####### Part 1
##Determining the gene expression in the individual samples using edgeR
## 1.Loading relevant modules
#install.packages("data.table")
library("data.table")
# install edgeR
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("edgeR")
#browseVignettes("edgeR")---to view the documentation
library("limma")  # limma must before edgeR
library("edgeR")

## 2. Set the current working directory
setwd("./all_htseq/")

## 3. Define variables
path_in ="htseq_all_count.csv"
group_in="sample_design.csv"
path_cpm ="cpm_normalised_htseq.csv"
path_voom ="voom_normalised_htseq.csv"

## 4. Read the counts matrix as a DGE object
x  = read.delim(path_in, header = T, row.names = 1, sep=",", stringsAsFactors = F)

## 5. Read the group data (sample design file)
sample_design = read.delim(group_in, header = T, sep=",", stringsAsFactors=F)
group = as.factor(c(sample_design$group))

## 6. Create a design object 
design = model.matrix(~ 0+group)
colnames(design) = c("normal","tumour")

## 7. modify the data in the sample_design file to match the data in x (i.e colnames have "." in one datasets and "-" in the other)
# the sample name of two files must be same, and same order

## 8. Replace all the "-" in the sting data with periods in group
sample_design$id = gsub("-", ".", sample_design$id)
rownames(sample_design) = sample_design$id
names(x) = sub("X*","",names(x))  # delete "X" in the coloum name of x (htseq_all_count.csv)

## 9. Create a dge list object
order = as.data.frame(t(sample_design))
x = x[names(order)]      #ordering colnames in X according to colnames in sample design
all(order(names(x))==order(names(order)))    #the answer should be TRUE
dge = DGEList(counts = x, group = group)     
dim(dge$counts)  # check the counts of each sample ----36924   139

## 10. Filter low-level expression genes
min = as.integer(dim(dge$counts)[2]*0.33)   
#keep samples that have at least an expression value that is greater than one
keep = rowSums(cpm(dge)>1) >= 1
dge  = dge[keep, , keep.lib.sizes=FALSE]  
dim(dge$counts)  # check the counts of each sample ----27632   139

## 11. Normalise the data
dge = calcNormFactors (dge)
#Plotting the data as an MDS object
#png("MDS.png")
plotMDS(dge, col = as.numeric(group))
#dev.off()

## 12. log cpm version
# way 1
logCPM = log2(cpm(dge)+1) 
#adding one before log conversion of cpm in order to avoid getting  infinite log0 values

# way 2-- edge R
logCPM_2 = cpm(dge, log= TRUE, prior.count=1)


##Writing output to a table
#write.csv(logCPM, file = path_cpm ,row.names = TRUE)

## Perform differential gene expression on logCPM vaues
fit = lmFit(logCPM, design)
fit = eBayes(fit, trend =TRUE)
results = topTable(fit, coef=ncol(design))

#Giving more weight to fold-changes in the gene ranking 
fit = lmFit(logCPM, design)
fit = treat(fit, lfc=log2(1.2) )
results_2 = topTreat(fit, coef=ncol(design))

## 13. voom version 
#Create a voom object
v = voom(dge, design, plot=TRUE)
#png("Voom.png")
#v = voom(dge, design, plot=TRUE)
#dev.off()
#Write output to a table
#write.csv(v$E, file = path_voom , row.names = TRUE)

####### Part 2
voomdataframe = log2((2^v$E) +1)
write.csv(voomdataframe, file = "voomdataframe.csv", row.names = TRUE)

library("ggplot2")
library("ggpubr")
library("tidyverse")
#library("rstatix")
library("grid")

setwd("all_htseq/")

### boxplot
data <- read.csv("voomdataframe-box (14 vs 125).csv",sep = ",")
stat.test <- compare_means(Expression ~ group,data=data,method = "t.test")
stat.test
data$group <- factor(data$group, levels = c("LNC14112-Normal","LNC14112-TNBC","REV3L-Normal","REV3L-TNBC"))
ggboxplot(data, x= "group", y ="Expression",fill="group") +
  stat_pvalue_manual(stat.test,y.position = 6.5,step.increase = 0.1,size=5,label = "p.signif") +
  scale_fill_manual(values = c("steelblue","yellow","steelblue","yellow")) +
  theme(legend.position = "bottom",plot.title = element_text(size = 20,hjust = 0.5),
        axis.title.x = element_blank()) +
  ggtitle("voom normalised data (14 normal vs. 125 diseased)")


####### Part 3
### correlation

data2 <- read.csv("voomdataframe-relationship (14 vs 125).csv",sep = ",")
sp <- ggplot(data=data2, aes(x=REV3L, y=LNC14112,colour=group)) + geom_point(size=3)
reg <- lm(LNC14112 ~ REV3L,data = data2)
coeff=coefficients(reg)
coeff
cor.test(data2$LNC14112,data2$REV3L,method = "pearson") # cor = 0.2588333  p-value = 0.002093
sp + geom_abline(intercept = 0.08373533,slope = 0.18019755, col="blue") + 
  ggtitle("voom normalised data (14 normal vs 125 diseased)") +
  theme(plot.title = element_text(size=25,hjust = 0.5),
        text = element_text(size = 15)) +
  annotation_custom(grobTree(textGrob("cor = 0.2588333  \np-value = 0.002093",x=0.1,y=0.95,hjust = 0,
                                      gp=gpar(col="black", fontsize=15, fontface="italic"))))


data <-  read.csv(file = "./correlation.csv",header = TRUE, sep = ",")  # file= correlation.csv

# data
shapiro.test(data$REV3L)        # p-value = 0.0001341

shapiro.test(data$l_LNC14112)   # p-value = 0.0003422
cor.test(data$l_LNC14112,data$REV3L,method = "pearson") # cor -0.430796 , p-value = 0.01232

line <- lm(l_LNC14112 ~ REV3L,data = data)
coeff = coefficients(line)
plot(data$REV3L,data$l_LNC14112,
     xlab = "REV3L", ylab = "l_LNC14112",cex.lab=1.5,cex.main=2.5,
     main="REV3L vs. l_LNC14112")
abline(line,col="blue")
legend(4.5,5.3,legend=c("cor = -0.430796","p-value = 0.01232"),
       box.lty = 0,cex=0.5)

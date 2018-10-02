####Instructions: Part A (of A and B)

##For this part of the lab you will use metagenomic and metatranscriptomic data generated from deep-sea sediment.
##The sediment was collected along a pollution gradient following the 2010 Deepwater Horizon oil spill.
##Distances along the seafloor increase away from the broken well head.
##Genomes were assembled, binned, and RNASeq reads were mapped to the genome database.
##RNASeq counts per coding DNA sequence (CDS) were obtained from the mapping file using HTSeq.
##You are provided with:
  #A table of RNASeq counts per CDS, 
  #Genome bin designations, coverage, G+C content, and completeness
  #Cursory CDS annotations
  #Sequencing library sizes

##Use this data to accomplish the following:
##1. Plot average and per bin coverages by sample (aka distance from oil spill)
##2. Plot contig data as contig coverage by %GC content
##3. Plot core gene data by genome coverage
##4. Calculate count data RPKM values
##5. Normalise RPKM data by coverage
##6. Determine differential gene expression between distal and near-well sites 
##7. Filter the dataset based for genes that are significantly differentially expressed (top 50)
##8. Sort the filtered dataset by log fold change
##9. Plot ordinations of all coverage unnormalised and normalised data
##10. Plot heatmaps of the top set of differentially expressed genes (using log transformed values and coverage un/normalised data)
##11. Plot heatmaps as above, but without unannotated proteins 


###The following protocol will help you complete the lab: parts A to D

### Part A.1: Packages and useful commands ###

##OPTIONAL: Check and re-set your working directory:
getwd()
setwd("/Users/kmhandley/Documents/1000-Courses-UoA/2018-BIONF-703/lab-2-703")
getwd()

##Set number of plots in panel:
par(mfrow=c(1,1))

##Get R version:
R.version.string
#[1] "R version 3.4.4 Patched (2018-03-19 r74624)"

###Tutorial on transcriptomic analyses

##DROPBOX LINK TO DATA FILES: https://www.dropbox.com/sh/ghknf2xmpxylgoi/AABjtKkz65TdkNlq7ABWm1EHa?dl=0
#r-unnorm-gulf-maxten.csv
#r-unnorm-gulf-maxten-annot.csv

##Check system packages and libraries:
sessionInfo()

##Install packages if not installed:
install.packages("gplots")
source("http://www.bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("EDASeq")

##Load (attach) libraries:
library("gplots")
library("edgeR")
require(EDASeq)



### Part A.2: data prep and exploration ###

##Load count data:
dataset <- read.csv("r-unnorm-gulf-maxten-annot.csv", row.names=1)
lib.size <- read.csv("r-library-size.csv", row.names=1)

##View first row of dataset:
dataset[1,]

##Select count columns, and view first row:
count.data <- dataset[,1:10]
count.data[1,]

##Normalize count data step 1: reads per kilobase of gene per million reads (RPKM)
#RPKM = (number of reads mapped to gene)*(1000/gene length)*(1000000/library size)
rpk <- (count.data)*(1000/dataset$gene.length); rpk[1:3,]; ncol(rpk)
colnames(rpk) <- gsub('count', 'rpk', colnames(rpk), fixed=TRUE); names(rpk)
rpkm <- rpk*(1000000/lib.size)[col(rpk)]; rpkm[1:3,]; ncol(rpkm)
colnames(rpkm) <- gsub('rpk', 'rpkm', colnames(rpkm), fixed=TRUE); names(rpkm)

##Normalize count data step 2: normalised counts per unit of genome coverage (e.g. thousand times)
#RPKM per thousand times coverage = RPKM * (1000/genome bin coverage)
coverage.data <- dataset[,c(11:20)]; coverage.data[1,]
rpkm.cov <- rpkm*(1000/coverage.data); rpkm.cov[1,]
colnames(rpkm.cov) <- gsub('rpkm', 'rpkm.cov', colnames(rpkm.cov), fixed=TRUE); names(rpkm.cov)

##Create and view list and objects in list of input data dataframes (list = collection of R objects):
datalist <- list()
datalist$lib.size <- lib.size
datalist$count.data <- count.data
datalist$rpk <- rpk
datalist$rpkm <- rpkm
datalist$rpkm.cov <- rpkm.cov
datalist; names(datalist)

##Plot coverage versus distance:
par(mfrow=c(1,1))
boxplot(log(coverage.data),margins=c(4,4))
##Plot coverage by bin:
par(mfrow=c(2,5))
sort <-dataset[order(dataset$bin.number),]
for (i in (11:20)) {boxplot(sort[,i]~sort$bin.number,ylim=c(0,450),main=colnames(sort[i]))}

##Plot bin coverage versus %GC content per gene:
par(mfrow=c(1,2))
plot(dataset$gc.gene,rowMeans(coverage.data)); abline(h=seq(20,80,by=20),v=seq(20,80,by=20),lty="dashed")
plot(dataset$core.gene.percent,rowMeans(coverage.data),pch=8,col="red"); abline(h=seq(20,80,by=20),v=seq(20,80,by=20),lty="dashed")



### Part A.3: Differential Gene Expression DGE analysis ###

##Group samples: (7 near-well; 3 distal)
sample.groups <- factor(c(1,1,1,1,1,1,1,2,2,2))

##Model matrix:
design <- model.matrix(~sample.groups)

##Create and view list and objects in list for differential gene expression (DGE):
list <- DGEList(counts=count.data,group=sample.groups)
list; names(list)

##Filter data by minimum number of counts: already done in this case (>=10 reads in >=1 samples)
#keep <- rowSums(list$counts>=10) >= 1
#table(keep)
#list <- list[keep, , keep.lib.sizes=FALSE]

##Change y$samples$lib.size to actual total library size, and view modified list: (not needed when offsets used)
nrow(datalist$lib.size)
list$samples$lib.size <- as.numeric(datalist$lib.size[1:10,])
list$samples

##RNASeq count data is relative. 
#Highly expressed genes can consume a substantial proportion of the total library size.
#This causes the remaining genes to be under-sampled in that sample.
#Therefore, you edgeR can apply scaling factors to minimise log-fold changes between samples.
#To do this use the calcNormFactors() function. 
list <- calcNormFactors(list) #Not applied when offsets include library size corrections.
list$samples

##Create offsets using EDASeq and inspect object:
#In EDASeq offset = log(ynorm + 0.1) − log(yraw + 0.1)
#where ynorm = "normalized counts" and yraw = "raw counts"
#EDASeq defines the offset as the log-ratio between "normalized and raw" counts. 
#edgeR expects the offset as the log-ratio between "raw and normalized" counts. 
#Therefore use -offst(offsetData) as the offset argument in edgeR.
list$offset <- offset(log(data.matrix(datalist$rpkm.cov)+0.1)-log(list$counts+0.1))
list$offset

##Estimate dispersions using edgeR:
#Overall (inherent) biological variability
#Must be accounted for to discriminate meaningful differences between conditions
disp <- estimateGLMCommonDisp(list$counts, design, offset=-offset(list$offset))

##Perform quasi-likelihood F-tests
#Use quasi-likelihood extension of a generalised linear model (negative binomial)
#The quasi negative binomial deals with over-dispersion - sequeezing estimates towards the mean
fit <- glmQLFit(list$counts, design, disp, offset=-offset(list$offset))
qlf <- glmQLFTest(fit, coef=2); names(qlf)

##Plot dispersions:
list <- estimateDisp(list, design); list$common.dispersion
#Plot gene-wise dispersions (with abundance trend):
par(mfrow=c(1,2)); plotBCV(list, main="Gene-wise dispersions")
#Plot QL dispersions:
plotQLDisp(fit, main="Raw and squeezed dispersions")

##Plot all the logFCs against average count size, highlighting the DE genes:
par(mfrow=c(1,1)); plotMD(qlf)
abline(h=c(-1,1), col="blue")
#Print number of differentially expressed genes:
summary(decideTests(qlf))

##Generate false discovery rate for dataset: 
#Testing for false positives arising from multiple comparisons
FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05)

##Combine output with dataset:
output<-cbind(dataset, rpk, rpkm, rpkm.cov, qlf$table, FDR)
names(output)


### Part A.4: Combine DGE info with alternative normalisation count data, etc ###

##Filter output to keep trustworthy, very differential expressed genes, and write output to file:
filt<-output[((output$logFC>=1.2 | output$logFC<=-1.2) & output$PValue<=0.05) & output$FDR<=0.05,]
filt[1,]; names(filt); nrow(output); nrow(filt)
write.table(filt, "r-filt-output.csv", sep=",")

##Filter output to keep trustworthy, very differential expressed genes, and excluding hypothetical annotations:
hyp.filt<-output[((output$logFC>=1.2 | output$logFC<=-1.2) & output$PValue<=0.05) & output$FDR<=0.05 & output$hypothetical=="known",]
hyp.filt[1,]; nrow(output); nrow(hyp.filt)
write.table(hyp.filt, "r-hyp-filt-output.csv", sep=",")

##Sort filtered data by log fold change, and select genes with the greatest change:
sort.filt <-filt[order(-(filt$logFC)),]
sortabs.filt <-filt[order(-abs(filt$logFC)),]
sortabs.filt$logFC; names(sortabs.filt); ncol(sortabs.filt)
top.sortabs.filt <- sortabs.filt[1:35,]
##Sort filtered data by p.value, and select genes with the greatest change:
sortpv.filt <-filt[order((filt$PValue)),]
sortpvabs.filt <-filt[order(abs(filt$PValue)),]
sortpvabs.filt$PValue; names(sortpvabs.filt); ncol(sortpvabs.filt)
top.sortpvabs.filt <- sortpvabs.filt[1:35,]

##Sort filtered data by log fold change, and select genes with the greatest change: annotated genes only
sort.hyp.filt <-hyp.filt[order(-(hyp.filt$logFC)),]
sortabs.hyp.filt <-hyp.filt[order(-abs(hyp.filt$logFC)),]
sortabs.hyp.filt$logFC; names(sortabs.hyp.filt); ncol(sortabs.hyp.filt)
top.sortabs.hyp.filt <- sortabs.hyp.filt[1:35,]

##Plot ordinations of filtered data sorted by log fold change
par(mfrow=c(1,2))
x <- as.matrix(sortabs.filt[,41:50])
plotMDS((log(x+1)), col=c(2,2,2,2,2,2,2,4,4,4),main="RPKM by log fold change")
x <- as.matrix(sortabs.filt[,51:60])
plotMDS((log(x+1)), col=c(2,2,2,2,2,2,2,4,4,4),main="RPKM.COV sorted by log fold change")

##Subset filtered data and plot as heatmaps: sorted by log fold change (annotated and hypothetical genes)
par(mfrow=c(1,1))
jet.colors <- colorRampPalette(c("white","blue")); pal <- jet.colors(256)
#All:
x <- as.matrix(sort.filt[,51:60])
hmap5 <- heatmap.2(log(x+1), col = pal, scale="none", margins=c(9,10), trace="none", density.info="none", keysize=1, symkey=FALSE, symbreaks=FALSE, main="Log RPKM.COV", labRow=FALSE)
#Top:
x <- as.matrix(top.sortabs.filt[,51:60])
hmap6 <- heatmap.2(log(x+1), col = pal, scale="none", margins=c(9,45), trace="none", density.info="none", keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="Log Top RPKM.COV")

##Subset filtered data and plot as heatmaps: sorted by log fold change (annotated genes)
#Top: Normalised
x <- as.matrix(top.sortabs.hyp.filt[,41:50])
hmap6 <- heatmap.2(log(x*100+1), col = pal, scale="none", margins=c(9,45), trace="none", density.info="none", keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="Log Top RPKM")
#Top: Normalised by coverage
x <- as.matrix(top.sortabs.hyp.filt[,51:60])
hmap6 <- heatmap.2(log(x+1), col = pal, scale="none", margins=c(9,45), trace="none", density.info="none", keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="Log Top RPKM.COV")



####Instructions: Part B (of A and B)

##For this part of the lab use metaproteomic data from a stream cyanobacterial bloom.
##Data represents biofilms were collected over 5 days (3 per day).
##Genomes were assembled, binned, and proteomic data were matched to the genome database.
##The data has undergone some filtering and normalisation steps already, but hasn't been corrected for genome coverage.
##Bin48 and bin38 are major and minor Phormidum (cyanobacterial) species in the biofilm.
##You are provided with:
##A table of coverage unnormalised data and genome coverage per sample per time point.


##Use this data to accomplish the following:
##1. Normalise to genome coverage
##2. Average the replicates per time point
##3. Filter data for photosynthesis related proteins
##4. Plot the coverage unnormalised and normalised photosynthesis data using heatmaps
##5. Do the same for nitrate-transport and P-transport

##Load count data:
dataset <- read.csv("omics-biofilm-unnorm.csv", row.names=1)
dataset[1,]

##Normalise by coverage:
cov <- as.matrix(dataset[,1:15])
unnorm <- as.matrix(dataset[,16:30])
norm <- unnorm*(10/cov); norm[1,]
colnames(norm) <- gsub('unnorm.', 'norm.', colnames(norm), fixed=TRUE)

##Average replicates:
#unnorm:
x <- as.matrix(unnorm); x[1:3,]
n <- 3; ave.unnorm <- apply(array(x, c(nrow(x), n, ncol(x)/n)), c(1, n), mean); ave.unnorm[1:n,]
colnames(ave.unnorm) <- c("ave.unnorm.day3","ave.unnorm.day6","ave.unnorm.day9","ave.unnorm.day12","ave.unnorm.day19")
ave.unnorm[1,]
#norm:
x <- as.matrix(norm); x[1:3,]
n <- 3; ave.norm <- apply(array(x, c(nrow(x), n, ncol(x)/n)), c(1, n), mean); ave.norm[1:n,]
colnames(ave.norm) <- c("ave.norm.day3","ave.norm.day6","ave.norm.day9","ave.norm.day12","ave.norm.day19")
ave.norm[1,]

##Combine datasets:
comb <- cbind(dataset,norm,ave.unnorm,ave.norm); comb[1,]
write.table(comb, "r-biofilm-normalised-output.csv", sep="\t")

##Filter output to keep photosystem proteins:
photo <-comb[comb$functional.category=="photosystem",]
photo[1,]; nrow(photo); ncol(photo)
write.table(photo, "r-biofilm-photosystems-output.csv", sep="\t")

##Heatmaps:
par(mfrow=c(1,1))
jet.colors <- colorRampPalette(c("white","blue")); pal <- jet.colors(256)
#Photosystems:
x <- as.matrix(photo[,58:62]); x[1,]
hmap1 <- heatmap.2(log(x+0.1), dendrogram="none", Rowv=FALSE, Colv=FALSE, col = pal, scale="row", margins=c(11,35), trace="none", density.info="none", keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="log unnormalised")
x <- as.matrix(photo[,53:57]); x[1,]
hmap2 <- heatmap.2(log(x+0.1), dendrogram="none", Rowv=FALSE, Colv=FALSE, col = pal, scale="row", margins=c(11,35), trace="none", density.info="none", keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="log normalised")

## Nitrate transport

##Filter output to keep P-transport proteins:
ptrans <-comb[comb$functional.category=="P-transport",]
ptrans[1,]; nrow(ptrans); ncol(ptrans)
write.table(ptrans, "r-biofilm-ptranssystems-output.csv", sep=",")

##Heatmaps:
par(mfrow=c(1,1))
jet.colors <- colorRampPalette(c("white","blue")); pal <- jet.colors(256)
#ptranssystems:
x <- as.matrix(ptrans[,58:62]); x[1,]
hmap1 <- heatmap.2(log(x+0.1), dendrogram="none", Rowv=FALSE, Colv=FALSE, col = pal, 
                   scale="row", margins=c(5,5), trace="none", density.info="none", 
                   keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="log unnormalised",
                   cexRow = 1, cexCol = 1)
x <- as.matrix(ptrans[,53:57]); x[1,]
hmap2 <- heatmap.2(log(x+0.1), dendrogram="none", Rowv=FALSE, Colv=FALSE, col = pal, 
                   scale="row", margins=c(11,11), trace="none", density.info="none", 
                   keysize=0.6, symkey=FALSE, symbreaks=FALSE, main="log normalised",
                   cexRow = 1, cexCol = 1)


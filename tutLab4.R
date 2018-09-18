library("vegan")
library("plyr")
library("RColorBrewer")
library("alphahull")
library("ggplot2")

# Read in the data
HPminus <- read.csv("HPminus.scaffold.coverage.csv", header = T)
HPplus <- read.csv("HPplus.scaffold.coverage.csv", header = T)
gc <- read.delim("assembly.gc.tab", header = T)
kmer <- read.delim("assembly.kmer.tab", header = T)
ess <- read.table("assembly.orfs.hmm.id.txt", header = F)
ess.tax <- read.delim("assembly.orfs.hmm.blast.tax.tab", header = F)
cons.tax <- read.delim("assembly.tax.consensus.txt", header = T)

colnames(kmer)[1] = "name"
colnames(ess) = c("name", "orf", "hmm.id")
colnames(ess.tax) = c("name", "orf", "phylum")
colnames(cons.tax) = c("name", "phylum", "tax.color", "all.assignments")

# Merge into single data frame d. Data on individual scaffolds
d <- as.data.frame(
  cbind(HPminus$Name, HPplus$Reference.length, gc$gc, HPminus$Average.coverage, 
        HPplus$Average.coverage), row.names = F)
colnames(d) = c("name", "length", "gc", "HPminus", "HPplus")    
d <- merge(d, cons.tax, by = "name", all = T)

# Tidying up the names
d$phylum <- sub("<phylum>", "", d$phylum)
d$phylum <- sub("unclassified Bacteria", "TM7", d$phylum)
d$phylum <- sub("/Chlorobi group", "", d$phylum)
d$phylum <- sub("Chlamydiae/", "", d$phylum)
d$phylum <- sub(" group", "", d$phylum)
d$phylum <- sub("Amoebozoa", NA, d$phylum)
d$phylum <- sub("Opisthokonta", NA, d$phylum)

# Essential gene dataframe
e <- merge(ess, d, by = "name", all.x = T)
e <- merge(e, ess.tax, by = c("name", "orf"), all.x = T)
e <- e[, -c(10, 11)]

# Stats on scaffolds
genome.stats <- matrix(NA, nrow = 0, ncol = 9)
colnames(genome.stats) <- c(
  "total.length", "# scaffolds", "mean.length", "max.length", "gc", 
  "HPminus", "HPplus", "tot.ess", "uni.ess")

calc.genome.stats <- function(x, y) 
  matrix(c(sum(x$length), nrow(x), round(mean(x$length), 1), max(x$length), 
           round(sum((x$gc * x$length))/sum(x$length), 1), 
           round(sum((x$HPminus * x$length))/sum(x$length), 1), 
           round(sum((x$HPplus * x$length))/sum(x$length), 1), nrow(y), 
           length(unique(y$hmm.id))), dimnames = list(colnames(genome.stats), 
                                                      ""))
extract <- function(x, a.def, v1, v2) {
  out <- {
  }
  for (i in 1:nrow(x)) {
    if (inahull(a.def, c(v1[i], v2[i]))) 
      out <- rbind(out, x[i, ])
  }
  return(out)
}

# Initial overview of the data:
calc.genome.stats(d, e)
# Subsets
ds <- subset(d, length > 5000)
es <- subset(e, length > 5000)

# Make a coverage plot
ggplot(ds, aes(x = HPminus, y = HPplus, color = gc, size = length)) + 
  scale_x_log10(limits=c(5,5000)) +
  scale_y_log10(limits=c(0.01,2000)) +
  xlab("Coverage (HP-)") +
  ylab("Coverage (HP+)") +
  geom_point(alpha = 0.5) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  scale_colour_gradientn(colours=c('red','green','blue'))

# colour by phylum
t <- 11
ds$tax.color[ is.na(ds$tax.color)] <- 0
for (i in 1:nrow(ds)){
  if (as.integer(ds$tax.color[i]) < t & as.integer(ds$tax.color[i]) > 0) {
    ds$tax.color[i] <- brewer.pal(8,'Paired')[as.integer(ds$tax.color[i])]
  } 
  else{
    ds$tax.color[i] <- NA
    ds$phylum[i] <- NA
  }
}

pcol<-cbind(unique(ds$tax.color)[-1],unique(ds$phylum)[-1])
pcol<-pcol[order(pcol[,2]),1]
ggplot(ds, aes(x = HPminus, y = HPplus, size = length, colour = phylum)) + 
  scale_x_log10(limits=c(5,5000)) +
  scale_y_log10(limits=c(0.01,2000)) +
  xlab("Coverage (HP-)") +
  ylab("Coverage (HP+)") +
  geom_point(alpha=0.1, colour = 'black') +
  geom_point(shape=1) +  
  scale_colour_manual(name="Phyla",values=pcol) +
  scale_size_area(name= "Scaffold length", max_size=20) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 5, shape = 19)))

# Genome extraction
x <- 'HPminus'
y <- 'HPplus'

plot(ds[,x], 
     ds[,y], 
     log="xy", 
     cex = sqrt(ds$length)/100, 
     pch=20, 
     col=rgb(0,0,0,0.1), 
     xlim = c(55,110),  
     ylim = c(0.5,10), 
     xlab = "Coverage HP-", 
     ylab = "Coverage HP+"
)

points(ds[,x], 
       ds[,y],
       cex = sqrt(ds$length)/100*0.7,
       col=ds$tax.color,
       lwd=2
)

def<-locator(100, type="p", pch=20)

#def<-{}
#def$x <- c(64,66,81,92,94,81,68,65)
#def$y <- c(2.0,6.6,7.7,3.9,1.4,1,1,1.4)

g1.selection.A <- ahull(def, alpha=100000)  

plot(g1.selection.A, col="black",add=T)

# Extract scaffolds and essential genes
g1.s.A<-extract(ds,g1.selection.A,ds[,x],ds[,y])
g1.e.A<-extract(es,g1.selection.A,es[,x],es[,y])

calc.genome.stats(g1.s.A, g1.e.A)

# Principal component analysis
rda <- rda(kmer[g1.s.A$name,2:ncol(kmer)],scale=T)
scores <- scores(rda,choices=1:5)$sites
# 
g1.s.B<-cbind(g1.s.A,scores)
g1.e.B<-merge(g1.e.A,g1.s.B[,c(1,9:13)],all.x=T,by="name")

rgb.c<- colorRampPalette(c('red','green','blue'))
rgb.a<-adjustcolor(rgb.c(max(d$gc)-min(d$gc)),alpha.f=0.2)
palette(rgb.a)

pairs(g1.s.B[,9:13], upper.panel=NULL, col = g1.s.B$gc-min(d$gc), 
      cex = sqrt(g1.s.B$length)/100, pch=20)

# Extract scaffolds using PCA instead
x <- 'PC1'
y <- 'PC2'

plot(g1.s.B[,x], 
     g1.s.B[,y], 
     cex = sqrt(g1.s.B$length)/100, 
     pch=20, 
     col=rgb(0,0,0,0.1), 
     xlab = x, 
     ylab = y
)

points(g1.s.B[,x], 
       g1.s.B[,y],
       cex = sqrt(g1.s.B$length)/100*0.7,
       col=g1.s.B$tax.color,
       lwd=1
)

def<-locator(100, type="p", pch=20)

#def<-{}
#def$x <- c(0.3740306,0.4839196,0.9084907,1.2431527,1.2781173,1.0733242,0.6537480,0.4689347,0.3690356)
#def$y <- c(0.28107380,1.31294166,1.94015545,1.99073721,1.33317436,0.39235367,0.04839772,0.02816501,0.22037569)

g1.selection.B <- ahull(def, alpha=100000)  

plot(g1.selection.B, col="black",add=T)

g1.s.C<-extract(g1.s.B,g1.selection.B,g1.s.B[,x],g1.s.B[,y])
g1.e.C<-extract(g1.e.B,g1.selection.B,g1.e.B[,x],g1.e.B[,y])

calc.genome.stats(g1.s.C, g1.e.C)
# Find duplicates
g1.d.C<-g1.e.C[which(duplicated(g1.e.C$hmm.id) | duplicated(g1.e.C$hmm.id, fromLast=TRUE)),] 
g1.d.C[order(g1.d.C$hmm.id),c(1,3,8)]

# Write bin to file
genome.stats<-rbind(genome.stats,t(calc.genome.stats(g1.s.C, g1.e.C)))
rownames(genome.stats)[nrow(genome.stats)]<-"genome 1"
show(genome.stats)
write.table(g1.s.C$name,file="genome1.txt",quote=F,row.names=F,col.names=F)


# bateriobates





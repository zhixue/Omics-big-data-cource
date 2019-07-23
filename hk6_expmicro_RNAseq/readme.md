<h2 style="text-align:center">Omics Big data

Homework 6  

</h2>

> 1. Using the microarray data you downloaded for the Question #3 in homework 2- 1, compare three different between-array normalization methods.
* keyword: hepatocellular carcinoma
* I use several between-array normalization methods, including **rma, mas5 ,li.wong and gcrma** here(but my data can not fit li.wong and gcrma due to same probe, Inf exists in mas5)
#### raw data
![E3BcNR.png](https://s2.ax1x.com/2019/04/29/E3BcNR.png)
![E3Bg41.png](https://s2.ax1x.com/2019/04/29/E3Bg41.png)
#### rma
![E3B7Ed.png](https://s2.ax1x.com/2019/04/29/E3B7Ed.png)
![E3BHUA.png](https://s2.ax1x.com/2019/04/29/E3BHUA.png)

#### Rcode
```R
source('http://bioconductor.org/biocLite.R')
library(GEOquery)
# get files
getGEOSuppFiles('GSE38199')
untar("GSE38199_RAW.tar",exdir = 'data_files')
cels <- list.files('data_files',pattern = '[gz]')
sapply(paste('data_files',cels,sep='/'),gunzip)

library(simpleaffy)
celfiles<-read.affy(covdesc = 'phenodata.txt', path = 'data_files')

# 1 plot
celfiles.rma<-rma(celfiles)
celfiles.mas5<-mas5(celfiles)

library(RColorBrewer)
cols<-rainbow(16)
boxplot(celfiles,col=cols)
hist(celfiles,col=cols)

library(affyPLM)
boxplot(celfiles.rma,col=cols)
hist(celfiles.rma,col=cols)

boxplot(celfiles.mas5,col=cols)
hist(celfiles.mas5,col=cols)
```

> 2. Based on the protocols of data analysis for microarray and RNA-Seq data, detect the DE genes (FDR<0.01) using your microarray and RNA-Seq data used in the Question 1. and then compare the top 100 significant DE genes obtained from microarray and RNA-Seq.

#### Rcode
```R
celfiles.filtered<-nsFilter(celfiles.rma,require.entrez = FALSE,remove.dupEntrez = FALSE)
filterEset<-exprs(celfiles.filtered$eset)

samples <- celfiles.rma$Target
samples <- as.factor(samples)
design<-model.matrix(~0 + samples)
colnames(design)<-c('pdgfc','wildtype')

library(limma)
fit <- lmFit(filterEset,design)
contrast.matrix<-makeContrasts(pdgfc_wildtype = pdgfc - wildtype,levels = design)
pdgfc_fits<-contrasts.fit(fit,contrast.matrix)

pdgfc_ebfit <-eBayes(pdgfc_fits)
head(topTable(pdgfc_ebfit,number = 100,coef = 1,p.value = 0.01))
```
* head of the results:
```
            logFC  AveExpr        t      P.Value    adj.P.Val        B
10548879 4.409561 7.334639 50.79324 4.075575e-20 7.083350e-16 34.24287
10480090 4.317457 7.229570 45.64248 2.510103e-19 2.181279e-15 32.94231
10352905 3.148716 6.340950 35.90773 1.464611e-17 8.484980e-14 29.73452
10356520 2.579317 7.048883 34.58412 2.763032e-17 1.200538e-13 29.20087
10470392 3.341377 7.408957 33.45460 4.840272e-17 1.682478e-13 28.72293
10594066 2.592172 6.022867 29.17562 4.850491e-16 1.139222e-12 26.69867)
```

> 3. Using your gene expression data above, calculate the correlation of gene expression measured by Microarray and RNA-Seq.

```R
exprset<-data.frame(exprs(celfiles.rma))
# correlation matrix
cor_expr<-cor(exprset)
# too big to show 
```

> 4. Using your RNA-Seq data above, calculate the correlation of gene expression measured by the raw count and RPKM/FPKM values.

* pipeline：
[![EYPWgU.md.png](https://s2.ax1x.com/2019/05/01/EYPWgU.md.png)](https://imgchr.com/i/EYPWgU)
### bash code
```bash
# mapping
tophat -p 8 -G ${genes}.gtf -o ${sample}_thout -–no-novel-juncs genome ${sample}_1.fq ${sample}_2.fq
# assemble expressed genes and transcripts
cufflinks -p 8 -o ${sample}_clout ${sample}_thout/accepted_hits.bam
# create a single merged transcriptome annotation
cuffmerge -g genes.gtf -s genome.fa -p 8 assemblies.txt
# Identify differentially expressed genes and transcript(control or treat)
cuffdiff -o diff_out -b genome.fa -p 8 –L C1,C2 -u merged_asm/merged.gtf \
./C1_R1_thout/accepted_hits.bam,./C1_R2_thout/accepted_hits.bam,./C1_R3_thout/
accepted_hits.bam \
./C2_R1_thout/accepted_hits.bam,./C2_R3_thout/accepted_hits.bam,./C2_R2_thout/
accepted_hits.bam
```

### Rcode

```R
# explore differential analysis results with cummerbund( in R)
library(cummeRbund)
cuff_data <- readCufflinks('diff_out')
csDensity(genes(cuff_data))
csScatter(genes(cuff_data), 'C1', 'C2')
csVolcano(genes(cuff_data), 'C1', 'C2')
Genes
mygene<- getGene(cuff_data,'regucalcin')
expressionBarplot(mygene)
expressionBarplot(isoforms(mygene))
# diff gene
gene_diff_data <- diffData(genes(cuff_data))
sig_gene_data <- subset(gene_diff_data, (significant == 'yes')) 
```
> 5. Using your microarray data, perform Hierarchical Clustering or K- means Clustering using the top 20% genes with most significant differential expression in the comparation between the case and the control groups. Show the heatmap andthefunctionalannotation(enrichmentanalysis)ofthemajorclusters. The codes and analysis workflow are also need to be provided.

![E3cC9A.png](https://s2.ax1x.com/2019/04/29/E3cC9A.png)
#### Rcode
```R
# 5
probeset.list<-topTable(pdgfc_ebfit,number = as.integer(0.2*length(pdgfc_ebfit$F)),coef = 1,p.value = 0.01)
selData<-filterEset[rownames(filterEset) %in% rownames(probeset.list),]

pdf(file='heatmap.pdf',labRow=c(''),col=topo.colors(16),cexCol=0.6)
heatmap(selData,labRow = c(''),col=topo.colors(16),cexCol = 0.6)
graphics.off()

#biocLite('hgu133plus2.db')
library(hgu133plus2.db)

genes.symbols <- getSYMBOL(rownames(probeset.list),'hgu133plus2')
results<-cbind(probeset.list,genes.symbols)
```
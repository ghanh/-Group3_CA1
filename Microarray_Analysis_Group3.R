# R version 3.6.2
# platform x86_64-w64-mingw32

## install required packages
# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery")
# BiocManager::install("limma")
# BiocManager::install("hgu133plus2.db")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("oligo")
# BiocManager::install("huex10sttranscriptcluster.db")

## load required packages
library(GEOquery)
library(affy)
library(limma)
library(hgu133plus2.db)
library(huex10sttranscriptcluster.db)
library(oligo)
library(RColorBrewer) #check the palettes here https://www.datanovia.com/en/blog/the-a-z-of-rcolorbrewer-palette/
library(dplyr)

## get data from depository
# download raw CEL files
# filePaths <- getGEOSuppFiles('GSE143150')

# set working directory to file that contains the .CEL.gz files
# setwd()
# getwd()

celFiles <- list.celfiles()

# parse GEO SOFT format file into an R data structure to make 
# access to each important part of GEO SOFT format easily accessible.
gse <- getGEO('GSE143150', GSEMatrix = FALSE)

# get_conditions takes in a GSM and gives the characteristics_ch1
get_conditions <- function(gsm) { 
      Meta(getGEO(gsm, GSEMatrix=FALSE))[['characteristics_ch1']] }

# a data frame is made using get_conditions function, contains useful metadata of each GSM of the GSE
pd <- data.frame((sapply(Meta(gse)[['sample_id']], get_conditions)), 
                 stringsAsFactors = TRUE)
pd <- t(pd)
colnames(pd) <- c('cell_line', 'genotype', 'treatment', 'age')
pd <- data.frame(pd, stringsAsFactors = TRUE)

# when renaming levels, need to be in lexicographical order!
levels(pd$cell_line) <- "HUVEC"
levels(pd$genotype) <- c('KD','WT')
levels(pd$treatment) <- c('DM','FM')
levels(pd$age) <- c('3d')
pd

### An AnnotatedDataFrame consists of two parts:
# (1) a collection of samples and the values of variables measured on those samples. 
# (2) a description of each variable measured. 
# Components can be accessed with pData() and varMetadata().
data <- read.celfiles(celFiles,phenoData=new("AnnotatedDataFrame",pd))

# background correction, normalisation and calculating expression with RMA
eset <- oligo::rma(data, background=TRUE, normalize=TRUE, target="core")

# plot expression data after RMA
affy::plotDensity(exprs(eset),xlab='log intensity',main="feature level densities after RMA",lwd=2)

## construct a linear model for the expression levels
# design and contrast matrices: 

cond <- as.factor(with(pData(eset),paste(genotype,treatment,sep=".")))
design_mat <- model.matrix(~ 0 + cond)
colnames(design_mat) <- levels(cond)

contrast <- makeContrasts(one = KD.DM - KD.FM,
                          two = WT.DM - WT.FM,
                          three = KD.DM - WT.DM,
                          four = KD.FM - WT.FM,
                          five = (KD.DM - WT.DM) - (KD.FM - WT.FM),
                          # , six = (KD.DM - WT.DM) + (KD.FM - WT.FM),
                          # seven = (WT.DM - WT.FM) +  (KD.DM - KD.FM),
                          # eight = WT.FM - KD.DM,
                          levels = design_mat
)

# fit the linear model
fit <- lmFit(eset, design_mat)
fit.contrast <- contrasts.fit(fit, contrast)
fit.bayes <- eBayes(fit.contrast)

# deg = differentially expressed genes
deg1 <- topTable(fit.bayes, coef = 'one', number = Inf)
deg2 <- topTable(fit.bayes, coef = 'two', number = Inf)
deg3 <- topTable(fit.bayes, coef = 'three', number = Inf)
deg4 <- topTable(fit.bayes, coef = 'four', number = Inf)
deg5 <- topTable(fit.bayes, coef = 'five', number = Inf)

## look at only interesting probesets for each contrast 
# by total number
# top_genes <- 50
# or by p-value and logfc
p = 0.05
logfc = 1.2
# degone <- topTable(fit.bayes, coef = 'one', number = Inf, p.value = p, lfc = logfc)
# degtwo <- topTable(fit.bayes, coef = 'two', number = Inf, p.value = p, lfc = logfc)
# degthree <- topTable(fit.bayes, coef = 'three', number = Inf, p.value = p, lfc = logfc)
# degfour <- topTable(fit.bayes, coef = 'four', number = Inf, p.value = p, lfc = logfc)
# degfive <- topTable(fit.bayes, coef = 'five', number = Inf, p.value = p, lfc = logfc)

comparison <- 'one'
interesting_genes1 <- topTable(fit.bayes, 
                              number = Inf, 
                              coef = comparison, 
                              p.value = p, lfc=logfc)
interesting_ids1 <- rownames(interesting_genes1)
eset_of_interest1 <- eset[interesting_ids1,]

# annotate gene names and database IDs for the interseting probesets
annot1 <- AnnotationDbi::select(huex10sttranscriptcluster.db, 
                               interesting_ids1, 
                               c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"),
                               keytype = "PROBEID"
)
rownames(interesting_genes1) <- annot1$SYMBOL
gene_list1 <- annot1$ENTREZID
gene_list1 <- gene_list1[!is.na(gene_list1)]

comparison <- 'two'
interesting_genes2 <- topTable(fit.bayes, 
                              number = Inf, 
                              coef = comparison, 
                              p.value = p, lfc=logfc)
interesting_ids2 <- rownames(interesting_genes2)
eset_of_interest2 <- eset[interesting_ids2,]
# annotate gene names and database IDs for the interseting probesets
annot2 <- AnnotationDbi::select(huex10sttranscriptcluster.db, 
                               interesting_ids2, 
                               c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME"),
                               keytype = "PROBEID"
)
rownames(interesting_genes2) <- annot2$SYMBOL
gene_list2 <- annot2$ENTREZID
gene_list2 <- gene_list2[!is.na(gene_list2)]



## make comparisons for EndMT markers between any two experimental groups
# retrieve marker gene names
marker_names <- c("TEK", "KDR", "PECAM1", "CDH5", "NOS3", "CNN1", "TAGLN", "VCAN", "SULF1", "CXCL1")
marker_IDs <- c("7010", "3791", "5175", "1003", "4846", "1264", "6876", "1462", "23213","2919")
panel <- AnnotationDbi::select(huex10sttranscriptcluster.db, 
                               marker_IDs, 
                               c("PROBEID", "GENENAME","SYMBOL"),
                               keytype = "ENTREZID"
)

# carry out t-tests on the expression levels of EndMT markers between any two groups
# and store the p-values
expression <- exprs(eset)
expression <- data.frame(expression, stringsAsFactors = TRUE)
dim(expression)
expression <- expression[rownames(expression)%in%panel$PROBEID,]
expression <- arrange(expression, match(rownames(expression), panel$PROBEID))
expression <- cbind(panel, expression)
WT.DM <- c("GSM4250986, GSM4250990", "GSM4250994")
KD.DM <- c("GSM4250987, GSM4250991", "GSM4250995")
WT.FM <- c("GSM4250988, GSM4250992", "GSM4250996")
KD.FM <- c("GSM4250989, GSM4250993", "GSM4250997")

WT.DMreps <- expression[, 4 + c(1,5,9)]
KD.DMreps <- expression[, 4 + c(2,6,10)]
WT.FMreps <- expression[, 4 + c(3,7,11)]
KD.FMreps <- expression[, 4 + c(4,8,12)]

WT.DM_FM <- vector()
for (i in 1:10) {
      WT.DM_FM <- c(WT.DM_FM, t.test(WT.DMreps[i,], WT.FMreps[i,])$p.value)
}
KD.DM_FM <- vector()
for (i in 1:10) {
      KD.DM_FM <- c(KD.DM_FM, t.test(KD.DMreps[i,], KD.FMreps[i,])$p.value)
}
DM.KD_WT <- vector()
for (i in 1:10) {
      DM.KD_WT <- c(DM.KD_WT, t.test(KD.DMreps[i,], WT.DMreps[i,])$p.value)
}
FM.KD_WT <- vector()
for (i in 1:10) {
      FM.KD_WT <- c(FM.KD_WT, t.test(KD.FMreps[i,], WT.FMreps[i,])$p.value)
}
# barplot for WT.DM vs WT.FM (two)
# get info on the probesets corresponding to the markers
panel_tbl <- deg2[rownames(deg2)%in%panel$PROBEID,]
panel_tbl
panel_tbl <- arrange(panel_tbl, match(rownames(panel_tbl), panel$PROBEID))
panel_tbl <- cbind(panel, panel_tbl)
panel_tbl

# construct a vector for the height of the bars
bar.vect <- vector()
for (i in panel_tbl$logFC) {
      bar.vect <- c(bar.vect, i)
}
names(bar.vect) <- panel_tbl$SYMBOL
bar <- barplot(bar.vect,
               cex.axis = 0.9, cex.names = 0.5,
               width = 0.25, border = NA,
               beside = TRUE,
               space = 0.3,
               legend.text = c("Endothelial", "Mesenchymal"),
               args.legend = list(bty = "n", x = "right", cex = 0.5, inset=c(-0.085,0), xpd = TRUE, fill = c("#32a18b", "#a13232")),
               xlim = c(0, 3.35), ylim = c(min(bar.vect) - 0.3, max(bar.vect) + 0.3),
               ylab = "log2 Fold Change",
               col = c(rep("#32a18b", 5), rep("#a13232", 6)),
               main = 'Effect of Differentiation Media on Wild-Type Samples', cex.main = 0.9
)

# make labels for the p-values and significance asterisks
sig <- vector()
for (i in WT.DM_FM) {
      if (i <= 0.05) {
            sig <- c(sig, 
                     paste0(paste0(replicate(floor(log10(1/i)), "*"), collapse = ""),
                            paste0("p=", signif(i, 3))))
      }
      else {
            sig <- c(sig, paste0("p=", signif(i, 3)))
      }
}

# ablinepiece() lets you draw a straight line of customisable length
ablinepiece <- function(a=NULL,b=NULL,reg=NULL,from,to,...){
      
      # Borrowed from abline
      if (!is.null(reg)) a <- reg
      
      if (!is.null(a) && is.list(a)) {
            temp <- as.vector(coefficients(a))
            
            if (length(temp) == 1) {
                  a <- 0
                  b <- temp
            }
            else {
                  a <- temp[1]
                  b <- temp[2]
            }
      }
      
      segments(x0=from,x1=to,
               y0=a+from*b,y1=a+to*b,...)
      
}

# draw a horizontal line through origin to represent 0 change level
ablinepiece(from = -1, to = 3.3, a = 0, b = 0)

# label p-value
text(bar, bar.vect + bar.vect/abs(bar.vect)*0.2, labels = as.character(sig), cex = 0.5)

### NOTE: the rest of the file only contains code similar to the chunk above to plot more bar charts ###
# barplot for KD.DM vs KD.FM (one)
panel_tbl <- deg1[rownames(deg1)%in%panel$PROBEID,]
panel_tbl
panel_tbl <- arrange(panel_tbl, match(rownames(panel_tbl), panel$PROBEID))
panel_tbl <- cbind(panel, panel_tbl)
panel_tbl

bar.vect <- vector()
for (i in panel_tbl$logFC) {
      bar.vect <- c(bar.vect, i)
}
names(bar.vect) <- panel_tbl$SYMBOL
bar <- barplot(bar.vect,
               cex.axis = 0.9, cex.names = 0.5,
               width = 0.25, border = NA,
               beside = TRUE,
               space = 0.3,
               legend.text = c("Endothelial", "Mesenchymal"),
               args.legend = list(bty = "n", x = "right", cex = 0.5, inset=c(-0.085,0), xpd = TRUE, fill = c("#32a18b", "#a13232")),
               xlim = c(0, 3.35), ylim = c(min(bar.vect) - 0.3, max(bar.vect) + 0.3),
               ylab = "log2 Fold Change",
               col = c(rep("#32a18b", 5), rep("#a13232", 6)),
               main = 'Effect of Differentiation Media on JMD2B-Silenced Samples', cex.main = 0.9
)
# make labels for the p-values and significance asterisks
sig <- vector()
for (i in KD.DM_FM) {
      if (i <= 0.05) {
            sig <- c(sig, 
                     paste0(paste0(replicate(floor(log10(1/i)), "*"), collapse = ""),
                            paste0("p=", signif(i, 3))))
      }
      else {
            sig <- c(sig, paste0("p=", signif(i, 3)))
      }
}

ablinepiece(from = -1, to = 3.3, a = 0, b = 0)
text(bar, bar.vect + bar.vect/abs(bar.vect)*0.2, labels = as.character(sig), cex = 0.5)


# barplot for KD.DM vs WT.DM (three)
panel_tbl <- deg3[rownames(deg3)%in%panel$PROBEID,]
panel_tbl
panel_tbl <- arrange(panel_tbl, match(rownames(panel_tbl), panel$PROBEID))
panel_tbl <- cbind(panel, panel_tbl)
panel_tbl

bar.vect <- vector()
for (i in panel_tbl$logFC) {
      bar.vect <- c(bar.vect, i)
}
names(bar.vect) <- panel_tbl$SYMBOL
bar <- barplot(bar.vect,
               cex.axis = 0.9, cex.names = 0.5,
               width = 0.25, border = NA,
               beside = TRUE,
               space = 0.3,
               legend.text = c("Endothelial", "Mesenchymal"),
               args.legend = list(bty = "n", x = "right", cex = 0.5, inset=c(-0.085,0), xpd = TRUE, fill = c("#32a18b", "#a13232")),
               xlim = c(0, 3.35), ylim = c(min(bar.vect) - 0.3, max(bar.vect) + 0.3),
               ylab = "log2 Fold Change",
               col = c(rep("#32a18b", 5), rep("#a13232", 6)),
               main = 'Effect of JMD2B Silencing in Differentiation Media', cex.main = 0.9
)
# make labels for the p-values and significance asterisks
sig <- vector()
for (i in DM.KD_WT) {
      if (i <= 0.05) {
            sig <- c(sig, 
                     paste0(paste0(replicate(floor(log10(1/i)), "*"), collapse = ""),
                            paste0("p=", signif(i, 3))))
      }
      else {
            sig <- c(sig, paste0("p=", signif(i, 3)))
      }
}

ablinepiece(from = -1, to = 3.3, a = 0, b = 0)
text(bar, bar.vect + bar.vect/abs(bar.vect)*0.1, labels = as.character(sig), cex = 0.5)


# barplot for KD.FM vs WT.FM (four)
panel_tbl <- deg4[rownames(deg4)%in%panel$PROBEID,]
panel_tbl
panel_tbl <- arrange(panel_tbl, match(rownames(panel_tbl), panel$PROBEID))
panel_tbl <- cbind(panel, panel_tbl)
panel_tbl

bar.vect <- vector()
for (i in panel_tbl$logFC) {
      bar.vect <- c(bar.vect, i)
}
names(bar.vect) <- panel_tbl$SYMBOL
bar <- barplot(bar.vect,
               cex.axis = 0.9, cex.names = 0.5,
               width = 0.25, border = NA,
               beside = TRUE,
               space = 0.3,
               legend.text = c("Endothelial", "Mesenchymal"),
               args.legend = list(bty = "n", x = "right", cex = 0.5, inset=c(-0.092,0), xpd = TRUE, fill = c("#32a18b", "#a13232")),
               xlim = c(0, 3.35), ylim = c(min(bar.vect) - 0.3, max(bar.vect) + 0.3),
               ylab = "log2 Fold Change",
               col = c(rep("#32a18b", 5), rep("#a13232", 6)),
               main = 'Effect of JMD2B Silencing in Full Media', cex.main = 0.9
)
# make labels for the p-values and significance asterisks
sig <- vector()
for (i in FM.KD_WT) {
      if (i <= 0.05) {
            sig <- c(sig, 
                     paste0(paste0(replicate(floor(log10(1/i)), "*"), collapse = ""),
                            paste0("p=", signif(i, 3))))
      }
      else {
            sig <- c(sig, paste0("p=", signif(i, 3)))
      }
}

ablinepiece(from = -1, to = 3.3, a = 0, b = 0)
text(bar, bar.vect + bar.vect/abs(bar.vect)*0.07, labels = as.character(sig), cex = 0.5)


# bar plot for (KD.DM - WT.DM) - (KD.FM - WT.FM) (five)
panel_tbl <- deg5[rownames(deg5)%in%panel$PROBEID,]
panel_tbl
panel_tbl <- arrange(panel_tbl, match(rownames(panel_tbl), panel$PROBEID))
panel_tbl <- cbind(panel, panel_tbl)
panel_tbl

bar.vect <- vector()
for (i in panel_tbl$logFC) {
      bar.vect <- c(bar.vect, i)
}
names(bar.vect) <- panel_tbl$SYMBOL
bar <- barplot(bar.vect,
               cex.axis = 0.9, cex.names = 0.5,
               width = 0.25, border = NA,
               beside = TRUE,
               space = 0.3,
               legend.text = c("Endothelial", "Mesenchymal"),
               args.legend = list(bty = "n", x = "right", cex = 0.5, inset=c(-0.1,0), xpd = TRUE, fill = c("#32a18b", "#a13232")),
               xlim = c(0, 3.35), ylim = c(min(bar.vect) - 0.2, max(bar.vect) + 0.2),
               ylab = "log2 Fold Change",
               col = c(rep("#32a18b", 5), rep("#a13232", 6)),
               main = 'Interaction between Type of Media and JMD2B Silencing', cex.main = 0.9
)
# make labels for the p-values and significance asterisks
sig <- vector()
for (i in p.adjust(panel_tbl$P.Value, method = "BH")) {
      if (i <= 0.05) {
            sig <- c(sig, 
                     paste0(paste0(replicate(floor(log10(1/i)), "*"), collapse = ""),
                            paste0("p=", signif(i, 3))))
      }
      else {
            sig <- c(sig, paste0("p=", signif(i, 3)))
      }
}

ablinepiece(from = -1, to = 3.3, a = 0, b = 0)

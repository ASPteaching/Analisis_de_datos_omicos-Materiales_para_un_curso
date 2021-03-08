## ----getCellTypesData, eval=FALSE-----------------------------------------------------------------------------------------
## library(affy)
## library(genefilter)
## affyPath<- "datos/celltypes/celfiles"
## cellTypesTargets<-  read.AnnotatedDataFrame("targets.txt", sep="\t",  path=affyPath)
## cellTypesTargets$filename = file.path(affyPath, row.names(cellTypesTargets))
## cellTypesRaw <- read.affybatch(cellTypesTargets$filename, phenoData=cellTypesTargets)
## eset_rma <- rma(cellTypesRaw)
## Filtered <- nsFilter(eset_rma)
## eset_rma_filtered <- Filtered[["eset"]]
## save(eset_rma, eset_rma_filtered, file="datos/celltypes/celltypes-normalized.rma.Rda")


## ----c7codec01, eval=T----------------------------------------------------------------------------------------------------
stopifnot(require(Biobase))
load (file="./datos/celltypes/celltypes-normalized.rma.Rda")
my.eset <- eset_rma_filtered
grupo_1 <- as.factor(pData(my.eset)$treat)
stopifnot(require(genefilter))
teststat <-rowttests(my.eset, "treat")
print(teststat[1:5,])


## ----sortBypvals, eval=T--------------------------------------------------------------------------------------------------
ranked <-teststat[order(teststat$p.value),]
print(ranked[1:5,])


## ----selectedGenes--------------------------------------------------------------------------------------------------------
selectedTeststat <- ranked[ranked$p.value < 0.01,]


## ----volcano0-------------------------------------------------------------------------------------------------------------
FC <- teststat$statistic
pVal<-teststat$p.value
Y <- -log(pVal)
plot(Y~FC)


## ---- eval=FALSE----------------------------------------------------------------------------------------------------------
## library(ggrepel)
## # plot adding up all layers we have seen so far
## volcanoP<- function (de,log2FoldChange,  pvalue,
##                      diffexpressed=NULL, col=NULL, delabel=NULL){
##   ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
##         geom_point() +
##         theme_minimal() +
##         geom_text_repel() +
##         scale_color_manual(values=c("blue", "black", "red")) +
##         geom_vline(xintercept=c(-0.6, 0.6), col="red") +
##         geom_hline(yintercept=-log10(0.05), col="red")
## }
## diffexpressed <- ifelse(abs(teststat$statistic)>10, TRUE, FALSE)
## label <- rep(NA, length(diffexpressed))
## label[diffexpressed] <- rownames(teststat)[diffexpressed]
## volcanoP (de=teststat, log2FoldChange=teststat$statistic,  pvalue=teststat$p.value,
##           diffexpressed = diffexpressed, delabel = label)


## ----adjustPvals----------------------------------------------------------------------------------------------------------
stopifnot(require(multtest))
procs <- c("Bonferroni","BH", "BY")
adjPvalues <- mt.rawp2adjp(teststat$p.value, procs)
names(adjPvalues)
ranked.adjusted<-cbind(ranked, adjPvalues$adjp)
head(ranked.adjusted)


## ----selectAdjusted-------------------------------------------------------------------------------------------------------
selectedAdjusted<-ranked.adjusted[ranked.adjusted$BY<0.001,]
nrow(selectedAdjusted)


## ---- echo=FALSE----------------------------------------------------------------------------------------------------------
stopifnot(require(Biobase))
load (file="./datos/celltypes/celltypes-normalized.rma.Rda")
my.eset <- eset_rma_filtered
targets <- pData(my.eset)
age.treat<- paste(targets$age,targets$treat, sep=".")
lev<-factor(age.treat, levels=unique(age.treat))
design <-model.matrix(~0+lev)
colnames(design)<-levels(lev)
rownames(design) <-rownames(targets)
print(design)


## ----setContrasts---------------------------------------------------------------------------------------------------------
require(limma)
cont.matrix <- makeContrasts (
      LPS.in.AGED=(Aged.LPS-Aged.MED),
      LPS.in.YOUNG=(Young.LPS-Young.MED),
      AGE=(Aged.MED-Young.MED),
      levels=design)
cont.matrix


## ----linearmodelfit,echo=F------------------------------------------------------------------------------------------------
require(limma)
fit<-lmFit(my.eset, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
save(fit.main, file="./datos/celltypes/celltypes-fit.main.Rda")


## ----print=FALSE, echo=TRUE-----------------------------------------------------------------------------------------------
topTab_LPS.in.AGED <- topTable (fit.main, number=nrow(fit.main), coef="LPS.in.AGED", adjust="fdr")
topTab_LPS.in.YOUNG <- topTable (fit.main, number=nrow(fit.main), coef="LPS.in.YOUNG", adjust="fdr")
topTab_AGE  <- topTable (fit.main, number=nrow(fit.main) , coef="AGE", adjust="fdr")


## ----volcano1-------------------------------------------------------------------------------------------------------------
coefnum = 1
opt <- par(cex.lab = 0.7)
volcanoplot(fit.main, coef=coefnum, highlight=10, names=fit.main$ID,
            main=paste("Differentially expressed genes",colnames(cont.matrix)[coefnum], sep="\n"))
abline(v=c(-1,1))
par(opt)


### R code from vignette source 'casoResuelto2.Rnw'

###################################################
### code chunk number 1: getPackages
###################################################
if (!(require("estrogen", character.only=T))){
    BiocManager::install("estrogen")
}
if (!(require("hgu95av2.db", character.only=T))){
  BiocManager::install("hgu95av2.db")
}



###################################################
### code chunk number 2: installBioC (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("estrogen")


###################################################
### code chunk number 3: preparaDirectorios
###################################################
workingDir <-getwd()
system("mkdir data")
system("mkdir results")
dataDir <-file.path(workingDir, "data")
resultsDir <- file.path(workingDir, "results")
setwd(workingDir)


###################################################
### code chunk number 4: fijaOpciones
###################################################
options(width=80)
options(digits=5)


###################################################
### code chunk number 5: estrogenDir
###################################################
require(estrogen)
estrogenDir <- system.file("extdata", package = "estrogen")
print(estrogenDir)


###################################################
### code chunk number 6: copyData (eval = FALSE)
###################################################
## system(paste ("cp ", estrogenDir,"/* ./data", sep=""))


###################################################
### code chunk number 7: affybatch.create
###################################################
require(Biobase)
require(affy)
sampleInfo <- read.AnnotatedDataFrame(file.path(estrogenDir,"targLimma.txt"), 
    header = TRUE, row.names = 1, sep="\t")
fileNames <- pData(sampleInfo)$FileName
rawData <- read.affybatch(filenames=file.path(estrogenDir,fileNames),
                          phenoData=sampleInfo)


###################################################
### code chunk number 8: wrongaffybatch.create
###################################################
require(affy)
setwd(estrogenDir)
rawData.wrong <- ReadAffy()
setwd(workingDir)


###################################################
### code chunk number 9: plotHist
###################################################
info <- data.frame(grupo=c(1,1,2,2,3,3,4,4))
sampleNames <- pData(rawData)$Target
hist(rawData, main="Signal distribution", col=info$grupo, lty=1:ncol(info))
legend (x="topright", legend=sampleNames , col=info$grupo, lty=1:ncol(info))


###################################################
### code chunk number 10: computeDeg
###################################################
deg<-AffyRNAdeg(rawData, log.it=T)
summaryAffyRNAdeg(deg) 
plotAffyRNAdeg(deg)
legend (x="bottomright", legend=sampleNames, col=1:nrow(info), lty=1:nrow(info), cex=0.7)


###################################################
### code chunk number 11: plotDeg (eval = FALSE)
###################################################
## plotAffyRNAdeg(deg)
## legend (x="bottomright", legend=sampleNames, col=1:nrow(info), lty=1:nrow(info), cex=0.7)


###################################################
### code chunk number 12: boxPlot
###################################################
### boxplot
boxplot(rawData, cex.axis=0.6, col=info$grupo, las=2, names=sampleNames)


###################################################
### code chunk number 13: plotDendro
###################################################
### La muestras del mismo grupo deber\'ian agruparse juntas
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plclust(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)
###


###################################################
### code chunk number 14: affyQCReport (eval = FALSE)
###################################################
## stopifnot(require(affyQCReport))
## QCReport(rawData,file=file.path(resultsDir,"QCReport.pdf"))


###################################################
### code chunk number 15: affyPLM
###################################################
stopifnot(require(affyPLM))
computePLM <- T
if(computePLM){
  Pset<- fitPLM(rawData)
  save(Pset, file=file.path(dataDir,"PLM.Rda"))
}else{
  load (file=file.path(dataDir,"PLM.Rda"))
}


###################################################
### code chunk number 16: plotRLE
###################################################
RLE(Pset, main = "Relative Log Expression", names=sampleNames, las=2, col=info$grupo+1, cex.axis=0.6,ylim=c(-5,5))


###################################################
### code chunk number 17: plotNUSE
###################################################
NUSE(Pset, main = "Normalized Unscaled Standard Errors", las=2, names=sampleNames, las=2, col=info$grupo+1, cex.axis=0.6, ylim=c(0.5,1.5))


###################################################
### code chunk number 18: cleanTheHouse
###################################################
rm(Pset)
gc()
detach("package:affyPLM")


###################################################
### code chunk number 19: normalization.rma
###################################################
stopifnot(require(affy))
normalize <- T
if(normalize){
  eset_rma <- rma(rawData)    
  save(eset_rma, file=file.path(dataDir,"normalized.Rda"))
}else{
  load (file=file.path(dataDir,"normalized.Rda"))
}


###################################################
### code chunk number 20: normBoxPlot
###################################################
boxplot(eset_rma,main="RMA", names=sampleNames, cex.axis=0.7, col=info$grupo+1,las=2)


###################################################
### code chunk number 21: compareNormalizations (eval = FALSE)
###################################################
## eset_mas5 <- mas5(rawData)  # Uses expresso (MAS 5.0 method) much slower than RMA!
## stopifnot(require(gcrma))
## eset_gcrma <- gcrma(rawData) # The 'library(gcrma)' needs to be loaded first.
## stopifnot(require(plier))
## eset_plier <- justPlier(rawData, normalize=T) # The 'library(plier)' needs to be loaded first.
## compara <-data.frame(RMA=exprs(eset_rma)[,1], MAS5 =exprs(eset_mas5)[,1],
##                     GCRMA=exprs(eset_gcrma)[,1], PLIER =exprs(eset_plier)[,1])
## pairs(compara)


###################################################
### code chunk number 22: filtratge
###################################################
require(genefilter)
require(hgu95av2.db)
filtrats <- nsFilter(eset_rma)
class(filtrats)
names(filtrats)
dim(exprs(filtrats$eset))


###################################################
### code chunk number 23: matDesign1
###################################################
design.1<-matrix(
c(1,1,0,0,0,0,0,0,
  0,0,1,1,0,0,0,0,
  0,0,0,0,1,1,0,0,
  0,0,0,0,0,0,1,1),
nrow=8,
byrow=F)
colnames(design.1)<-c("neg10h", "est10h", "neg48h", "est48h")
rownames(design.1) <-  c("low10A", "low10B", "hi10A" , "hi10B",  "low48A", "low48B", "hi48A" , "hi48B") 
print(design.1)


###################################################
### code chunk number 24: selectLimma
###################################################
require(Biobase)
if (!exists("eset_rma"))  load(file.path(dataDir, "normalized.rda"))
targets <- pData(eset_rma)
stopifnot(require(limma))
lev<-factor(targets$Target, levels=unique(targets$Target))
design <-model.matrix(~0+lev)
colnames(design)<-levels(lev)
rownames(design) <-rownames(targets)
print(design)


###################################################
### code chunk number 25: setContrasts
###################################################
cont.matrix <- makeContrasts (
      Estro10=(est10h-neg10h),
      Estro48=(est48h-neg48h),
      Tiempo=(neg48h-neg10h),
      levels=design)
cont.matrix


###################################################
### code chunk number 26: linearmodelfit
###################################################
require(limma)
fit<-lmFit(eset_rma, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)


###################################################
### code chunk number 27: casoResuelto2.Rnw:511-514
###################################################
topTabEstro10 <- topTable (fit.main, number=nrow(fit.main), coef="Estro10", adjust="fdr")
topTabEstro48 <- topTable (fit.main, number=nrow(fit.main), coef="Estro48", adjust="fdr")
topTabTiempo  <- topTable (fit.main, number=nrow(fit.main) , coef="Tiempo", adjust="fdr")


###################################################
### code chunk number 28: volcano1
###################################################
coefnum = 1
opt <- par(cex.lab = 0.7)
volcanoplot(fit.main, coef=coefnum, highlight=10, names=rownames(fit.main), 
            main=paste("Differentially expressed genes",colnames(cont.matrix)[coefnum], sep="\n"))
abline(v=c(-1,1))
par(opt)


###################################################
### code chunk number 29: anota1
###################################################
require(hgu95av2.db)
hgu95av2()


###################################################
### code chunk number 30: genesEstro
###################################################
top5 <-rownames(topTabEstro10)[1:5]
cat("Usando mget\n")
geneSymbol5.1 <- unlist(mget(top5, hgu95av2SYMBOL))
geneSymbol5.1
cat("Usando toTable\n")
genesTable<- toTable(hgu95av2SYMBOL)
rownames(genesTable) <-  genesTable$probe_id
genesTable[top5, 2]
cat("Usando getSYMBOL\n")
require(annotate)
geneSymbol5.3 <- getSYMBOL(top5, "hgu95av2.db")
geneSymbol5.3


###################################################
### code chunk number 31: decideTests.1
###################################################
stopifnot(require(annotate))
anotPackage <- annotation(eset_rma)
fit.Symbols <- getSYMBOL (rownames(fit.main), anotPackage)
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.01)


###################################################
### code chunk number 32: resumeDecideTests
###################################################
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,]
print(summary(res))


###################################################
### code chunk number 33: venn1
###################################################
vennDiagram (res.selected[,1:3], main="Genes in common #1", cex=0.9)


###################################################
### code chunk number 34: prepareData
###################################################
probeNames<-rownames(res)
probeNames.selected<-probeNames[sum.res.rows!=0]
exprs2cluster <-exprs(eset_rma)[probeNames.selected,]


###################################################
### code chunk number 35: plotHeatMap1
###################################################
color.map <- function(horas) { if (horas< 20) "yellow" else "red" }
grupColors <- unlist(lapply(pData(eset_rma)$time.h, color.map))
heatmap(exprs2cluster, col=rainbow(100), ColSideColors=grupColors, cexCol=0.9)


###################################################
### code chunk number 36: plotHeatMap2
###################################################
color.map <- function(horas) { if (horas< 20) "yellow" else "red" }
grupColors <- unlist(lapply(pData(eset_rma)$time.h, color.map))

require("gplots")
heatmap.2(exprs2cluster, 
          col=bluered(75), scale="row",
          ColSideColors=grupColors, key=TRUE, symkey=FALSE, 
          density.info="none", trace="none", cexCol=1)



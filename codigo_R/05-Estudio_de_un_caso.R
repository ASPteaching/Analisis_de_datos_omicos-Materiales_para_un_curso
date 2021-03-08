## ----getPackages----------------------------------------------------------------------------------------------------------
if (!(require("estrogen", character.only=T))){
    BiocManager::install("estrogen")
    }


## ----installBioC, eval=F--------------------------------------------------------------------------------------------------
## BiocManager::install()


## ----preparaDirectorios---------------------------------------------------------------------------------------------------
workingDir <-getwd()
if (!file.exists("datos")) system("mkdir datos")
if (!file.exists("datos/estrogen")) system("mkdir datos/estrogen")
if (!file.exists("results")) system("mkdir results")
dataDir <-file.path(workingDir, "datos/estrogen")
resultsDir <- file.path(workingDir, "results")
setwd(workingDir)


## ----fijaOpciones---------------------------------------------------------------------------------------------------------
options(width=80)
options(digits=5)


## ----estrogenDir----------------------------------------------------------------------------------------------------------
library(estrogen)
estrogenDir <- system.file("extdata", package = "estrogen")
print(estrogenDir)


## ----copyData, eval=F-----------------------------------------------------------------------------------------------------
## system(paste ("cp ", estrogenDir,"/* ./data", sep=""))


## ----affybatch.create-----------------------------------------------------------------------------------------------------
library(Biobase)
library(affy)
sampleInfo <- read.AnnotatedDataFrame(file.path(estrogenDir,"targLimma.txt"),
    header = TRUE, row.names = 1, sep="\t")
fileNames <- pData(sampleInfo)$FileName
rawData <- read.affybatch(filenames=file.path(estrogenDir,fileNames),
                          phenoData=sampleInfo)


## ----wrongaffybatch.create, eval=T----------------------------------------------------------------------------------------
library(affy)
sampleInfoWithBad <- read.AnnotatedDataFrame(file.path(dataDir,"phenoDataWithBad.txt"),
    header = TRUE, row.names = NULL, sep="\t")
fileNames <- pData(sampleInfoWithBad)$FileName
rawData.wrong <- read.affybatch(filenames=file.path(dataDir,fileNames),
                          phenoData=sampleInfoWithBad)


## ----plotHist,echo=F,fig=T------------------------------------------------------------------------------------------------
info <- data.frame(grupo=c(1,1,2,2,3,3,4,4))
sampleNames <- pData(rawData)$Target
hist(rawData, main="Signal distribution", col=info$grupo, lty=1:ncol(info))
legend (x="topright", legend=sampleNames , col=info$grupo, lty=1:ncol(info))


## ----computeDeg, echo=F---------------------------------------------------------------------------------------------------
deg<-AffyRNAdeg(rawData, log.it=T)
summaryAffyRNAdeg(deg)
# plotAffyRNAdeg(deg)
# legend (x="bottomright", legend=sampleNames, col=1:nrow(info), lty=1:nrow(info), cex=0.7)


## ----plotDeg,echo=F, eval=F-----------------------------------------------------------------------------------------------
## plotAffyRNAdeg(deg)
## legend (x="bottomright", legend=sampleNames, col=1:nrow(info), lty=1:nrow(info), cex=0.7)


## ----boxPlot,echo=F,fig=T-------------------------------------------------------------------------------------------------
boxplot(rawData, cex.axis=0.6, col=info$grupo, las=2, names=sampleNames)


## ----plotDendro2, echo=F,fig=T--------------------------------------------------------------------------------------------
### La muestras del mismo grupo deberían agruparse juntas
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)
###


## ----arrayQualityMetrics, eval=F------------------------------------------------------------------------------------------
## if(!require(arrayQualityMetrics)) BiocManager::install("arrayQualityMetrics")
## library(arrayQualityMetrics)
## arrayQualityMetrics(rawData, outdir = "Informe_de_calidad_para_los_datos_del_caso_ESTROGENO")
## arrayQualityMetrics(rawDataBad, outdir = "Informe_de_calidad_(2)_para_los_datos_del_caso_ESTROGENO")


## ----normalizationrma-----------------------------------------------------------------------------------------------------
library(affy)
eset_rma <- rma(rawData)


## ----normalizationrma2, eval=FALSE----------------------------------------------------------------------------------------
## library(affy)
## normalize <- T
## if(normalize){
##   eset_rma <- rma(rawData)
##   save(eset_rma, file=file.path(dataDir,"estrogen-normalized.Rda"))
## }else{
##   load (file=file.path(dataDir,"estrogen-normalized.Rda"))
## }


## ----normBoxPlot, fig=TRUE------------------------------------------------------------------------------------------------
boxplot(eset_rma,main="Boxplot for RMA-normalized expression values ",
        names=sampleNames, cex.axis=0.7, col=info$grupo+1,las=2)


## ----compareNormalizations, echo=T, eval=F--------------------------------------------------------------------------------
## eset_mas5 <- mas5(rawData)  # Uses expresso (MAS 5.0 method) much slower than RMA!
## stopifnot(require(gcrma))
## eset_gcrma <- gcrma(rawData) # The 'library(gcrma)' needs to be loaded first.
## stopifnot(require(plier))
## eset_plier <- justPlier(rawData, normalize=T) # The 'library(plier)' needs to be loaded first.
## compara <-data.frame(RMA=exprs(eset_rma)[,1], MAS5 =exprs(eset_mas5)[,1],
##                     GCRMA=exprs(eset_gcrma)[,1], PLIER =exprs(eset_plier)[,1])
## pairs(compara)


## ----filtraje-------------------------------------------------------------------------------------------------------------
require(genefilter)
filtered <- nsFilter(eset_rma, require.entrez=TRUE,
         remove.dupEntrez=TRUE, var.func=IQR,
         var.cutoff=0.5, var.filter=TRUE,
         filterByQuantile=TRUE, feature.exclude="^AFFX")


## ----filtrado-------------------------------------------------------------------------------------------------------------
names(filtered)
class(filtered$eset)
print(filtered$filter.log)
eset_filtered <-filtered$eset


## ----saveData-------------------------------------------------------------------------------------------------------------
save(eset_rma, eset_filtered, file=file.path(resultsDir, "estrogen-normalized.Rda"))


## ----matDesign1-----------------------------------------------------------------------------------------------------------
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


## ----selectLimma, echo=F--------------------------------------------------------------------------------------------------
require(Biobase)
if (!exists("eset_rma"))  load(file.path(dataDir, "estrogen-normalized.Rda"))
targets <- pData(eset_rma)
stopifnot(require(limma))
lev<-factor(targets$Target, levels=unique(targets$Target))
design <-model.matrix(~0+lev)
colnames(design)<-levels(lev)
rownames(design) <-rownames(targets)
print(design)


## ----setContrasts2--------------------------------------------------------------------------------------------------------
cont.matrix <- makeContrasts (
      Estro10=(est10h-neg10h),
      Estro48=(est48h-neg48h),
      Tiempo=(neg48h-neg10h),
      levels=design)
cont.matrix


## ----linearmodelfit2------------------------------------------------------------------------------------------------------
require(limma)
fit<-lmFit(eset_rma, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)


## ----print=FALSE, echo=TRUE-----------------------------------------------------------------------------------------------
topTabEstro10 <- topTable (fit.main, number=nrow(fit.main), coef="Estro10", adjust="fdr")
topTabEstro48 <- topTable (fit.main, number=nrow(fit.main), coef="Estro48", adjust="fdr")
topTabTiempo  <- topTable (fit.main, number=nrow(fit.main) , coef="Tiempo", adjust="fdr")


## ----volcano2,echo=F, fig=T-----------------------------------------------------------------------------------------------
coefnum = 1
opt <- par(cex.lab = 0.7)
volcanoplot(fit.main, coef=coefnum, highlight=10, names=fit.main$ID,
            main=paste("Differentially expressed genes",colnames(cont.matrix)[coefnum], sep="\n"))
abline(v=c(-1,1))
par(opt)


## ----decideTests.1, echo=F------------------------------------------------------------------------------------------------
stopifnot(require(annotate))
anotPackage <- annotation(eset_rma)
fit.Symbols <- getSYMBOL (rownames(fit.main), anotPackage)
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.01)


## ----resumeDecideTests----------------------------------------------------------------------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,]
print(summary(res))


## ----venn1,fig=T----------------------------------------------------------------------------------------------------------
vennDiagram (res.selected[,1:3], main="Genes in common #1", cex=0.9)


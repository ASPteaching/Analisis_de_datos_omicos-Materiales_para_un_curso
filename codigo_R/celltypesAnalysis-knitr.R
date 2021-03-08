## ----include=FALSE-------------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)


## ----loadPackages--------------------------------------------------------------
installPackages=TRUE
if(installPackages){
  source("http://bioconductor.org/biocLite.R")
  if (!(require(Biobase))) biocLite("Biobase")
  if (!(require(affy))) biocLite("affy")
  if (!(require(arrayQualityMetrics))) biocLite("arrayQualityMetrics")
  if (!(require("mouse4302.db"))) biocLite("mouse4302.db")
  if (!(require("multtest"))) biocLite("multtest")
  if (!(require("limma"))) biocLite("limma")
  if (!(require("GOstats"))) biocLite("GOstats")
}else{
  require(Biobase) 
  require(affy)
  require(arrayQualityMetrics)
  require("mouse4302.db") 
  require("multtest")
  require("limma")
  require("GOstats")
}


## ----celltypesCreate-----------------------------------------------------------
workingDir<-getwd()
dataDir <-file.path(workingDir, 'data')
resultsDir <- file.path(workingDir,'results' )
require(Biobase)
require(affy)
sampleInfo <- read.AnnotatedDataFrame(file.path(dataDir,"targets.txt"), 
    header = TRUE, row.names = 1, sep="\t")
celltypes.experimentInfo <- new("MIAME", name="LPS_Experiment",
          lab="National Cancer Institute",
          contact="Lakshman Chelvaraja",
          title="Molecular basis of age associated cytokine dysregulation in LPS stimulated macrophages ",
          url="http://www.jleukbio.org/cgi/content/abstract/79/6/1314")
info <- pData(sampleInfo)
fileNames <-rownames(info)
sampleNames <- unlist(strsplit(fileNames,"\\."))[seq(1,23, by=2)]
rawData <- read.affybatch(filenames=file.path(dataDir,fileNames),
                          phenoData=sampleInfo, 
                          description=celltypes.experimentInfo)
save(info, sampleInfo, rawData, file=file.path(dataDir,"rawData.Rda"))


## ----plotHist,echo=F-----------------------------------------------------------
hist(rawData, main="Signal distribution", col=info$grupo, lty=1:ncol(info))
legend (x="topright", legend=sampleNames, col=info$grupo, lty=1:ncol(info))


## ----plotDeg,echo=F------------------------------------------------------------
deg<-AffyRNAdeg(rawData, log.it=T)
summaryAffyRNAdeg(deg) 
plotAffyRNAdeg(deg,col=1:nrow(info))
legend (x="bottomright", legend=sampleNames, col=1:nrow(info), lty=1:nrow(info), cex=0.7)


## ----plotBox,echo=F------------------------------------------------------------
### boxplot
boxplot(rawData, cex.axis=0.6, col=info$grupo, las=2, names=sampleNames)


## ----plotDendro, echo=F--------------------------------------------------------
### La muestras del mismo grupo deber?an agruparse juntas
clust.euclid.average <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of samples",  hang=-1)
###


## ----plotPCA, echo=F-----------------------------------------------------------
## ----funPCA--------------------------------------------------------------
#Definici?n de la funci?n de PCA-2d
#vigilar los "ylim". Habr?a que definirlos en cada caso
plotPCA <- function ( X, labels=NULL, colors=NULL, var = "", dataDesc="", 
                     scale=FALSE, formapunts=NULL, myCex=NULL,...)
{
  pcX<-prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
 if (is.null(colors)) colors=colores
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab,
       #xlim=c(min(pcX$x[,1])-5,max(pcX$x[,1])+5),
       #ylim=c(-35,50),
       pch=formapunts, col=colors)
        text(pcX$x[,1],pcX$x[,2],labels,pos=3,cex=myCex, col=colors)
  title(dataDesc, cex=0.2)
}
plotPCA(X=exprs(rawData), labels=sampleNames, dataDesc="Principal Components (Raw Data)", 
        var = "info$grupo", myCex=0.6, colors = info$grupo, formapunts=info$grupo)


## ----arrayQM1, eval=FALSE------------------------------------------------------
## require(arrayQualityMetrics)
## arrayQualitymetrics(rawData)


## ----readData3-----------------------------------------------------------------
stopifnot(require(affy))
if (!exists("rawData")) load(file=file.path (dataDir, "rawData.Rda"))


## ----normalization.rma---------------------------------------------------------
stopifnot(require(affy))
normalize <- T
if(normalize){
  eset_rma <- rma(rawData)    # Creates expression values using RMA method. 
  save(eset_rma, file=file.path(dataDir,"normalized.Rda"))
}else{
  load (file=file.path(dataDir,"normalized.rma.Rda"))
}


## ----normBoxPlot---------------------------------------------------------------
boxplot(eset_rma,main="RMA", names=sampleNames, cex.axis=0.7, col=info$grupo+1,las=2)


## ----compareNormalizations, echo=T, eval=F-------------------------------------
## eset_mas5 <- mas5(rawData)  # Uses expresso (MAS 5.0 method) much slower than RMA!
## stopifnot(require(gcrma))
## eset_gcrma <- gcrma(rawData) # The 'library(gcrma)' needs to be loaded first.
## stopifnot(require(plier))
## eset_plier <- justPlier(rawData, normalize=T) # The 'library(plier)' needs to be loaded first.
## compara <-data.frame(RMA=exprs(eset_rma)[,1], MAS5 =exprs(eset_mas5)[,1],
##                     GCRMA=exprs(eset_gcrma)[,1], PLIER =exprs(eset_plier)[,1])
## pairs(compara)


## ----relativeImportance, echo=F------------------------------------------------
plotPCA(X=exprs(eset_rma), labels=sampleNames, dataDesc="Principal Components (Normalized data)", 
        var = "info$grupo", myCex=0.6, colors = info$grupo, formapunts=info$grupo)


## ----filter1-------------------------------------------------------------------
library(genefilter)
log_intensity_threshold <-10
numsamples <- nrow(pData(eset_rma))/2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)
which1 <- genefilter(exprs(eset_rma), ffun1)
sum(which1)


## ----filter2-------------------------------------------------------------------
 sds<- apply(exprs(eset_rma), 1, sd)
 variability_threshold <- quantile (sds,0.9)
 f2 <-function (x) if (sd(x)< variability_threshold) return(FALSE) else return(TRUE)
 ffun2<-filterfun(f2)
 which2 <- genefilter(exprs(eset_rma), ffun2)
 sum (which2)


## ----commongenes---------------------------------------------------------------
probes <-rownames(exprs(eset_rma))
length(intersect(probes[which1],probes[which2]))


## ----filtering3----------------------------------------------------------------
require(genefilter)
require(mouse4302.db)
annotation(eset_rma)<-"mouse4302.db"
filtered <- nsFilter(eset_rma, require.entrez=TRUE,
         remove.dupEntrez=TRUE, var.func=IQR,
         var.cutoff=0.5, var.filter=TRUE,
         filterByQuantile=TRUE, feature.exclude="^AFFX")


## ----filtrado------------------------------------------------------------------
names(filtered)
class(filtered$eset)
print(filtered$filter.log)
eset_filtered <-filtered$eset
dim(exprs(eset_filtered))


## ----storeData,echo=F, results='hide'------------------------------------------
save(eset_rma, eset_filtered,
      file=file.path(dataDir, "normalized.filtered.Rda"))


## ----loadAgedData--------------------------------------------------------------
stopifnot(require(Biobase)) #library(Biobase)
load(file=file.path(dataDir, "normalized.filtered.Rda"))
my.eset <- eset_filtered[,pData(eset_filtered)$age=="Aged"]


## ----teststat------------------------------------------------------------------
stopifnot(require(genefilter))
# we need a class vector made of 0 qand 1s
teststat <-rowttests(my.eset, "treat")


## ----showTeststat, echo=F------------------------------------------------------
print(teststat[1:10,])


## ----hist2---------------------------------------------------------------------
hist(teststat$statistic)


## ----qqnorm1-------------------------------------------------------------------
qqnorm(teststat$statistic) 
qqline(teststat$statistic)


## ----sortStatistics------------------------------------------------------------
anotPackage="mouse4302.db"
require(anotPackage, character.only = T)
require(annotate)
probenames <- rownames(exprs(my.eset))
geneNames <-unlist(getSYMBOL(probenames, anotPackage))
topDown<-order(teststat$p.value)
ranked<-data.frame(gene=geneNames[topDown], t=teststat$statistic[topDown], 
  foldChg =teststat$dm[topDown] , pvalue=teststat$p.value[topDown])
print(top25<-ranked[1:10,])


## ----volcanoPlot, echo=F-------------------------------------------------------
x <-ranked$foldChg
y <- -log(ranked$pvalue)
plot(x, y, xlab="Fold Change", ylab ="-log(p-value", main="Volcano Plot")
abline(v=-2);abline(v=2);
text (x[1:25], y[1:25],ranked$gene[1:25],cex=0.7)


## ----selectNaif----------------------------------------------------------------
selectedNaif <-ranked[ranked$pvalue<0.001,]
nrow(selectedNaif)


## ----adjustPvalues-------------------------------------------------------------
stopifnot(require(multtest))
#procs <- c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH", "BY")
procs <- c("Bonferroni","BH", "BY")
adjPvalues <- mt.rawp2adjp(ranked$pvalue, procs)
names(adjPvalues)
ranked.adjusted<-cbind(ranked[,c(1,4)], adjPvalues$adjp[,-1])
ranked.adjusted[1:10,]


## ----selectedAdjusted----------------------------------------------------------
selectedAdjusted<-ranked.adjusted[ranked.adjusted$BY<0.01,]
nrow(selectedAdjusted)


## ----selectLimma, echo=F-------------------------------------------------------
stopifnot(require(limma))
lev<-as.factor(pData(my.eset)$treat)
design <-model.matrix(~0+lev)
colnames(design)<-paste("Grp",as.character(unique(lev)),sep="")
rownames(design) <-rownames(pData(my.eset))
print(design)


## ----setContrasts--------------------------------------------------------------
cont.matrix <- makeContrasts (LPSvsMED=(GrpLPS-GrpMED), levels=design)
cont.matrix


## ----linearmodelfit,echo=F-----------------------------------------------------
fit<-lmFit(my.eset, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)


## ----topTabs1------------------------------------------------------------------
topTab<-topTable(fit.main, number=nrow(fit.main), adjust="fdr")
geneNames <-getSYMBOL(rownames(topTab), anotPackage)
topTab<-cbind(gene=geneNames, topTab)
topTab.selected<-topTab[topTab$B>10,]
print(topTab.selected[1:10,])
write.table(topTab.selected$ID, file="ProbesList.txt")
write.table(topTab.selected, file="topTab.selected.txt", sep="\t")


## ----volcano2------------------------------------------------------------------
coefnum = 1
opt <- par(cex.lab = 0.7)
symbols.main<-geneNames[unlist(fit.main$genes)]
volcanoplot(fit.main, coef=coefnum, highlight=10, names=symbols.main)
abline(v=c(-1,1))
par(opt)


## ----commonGenes---------------------------------------------------------------
c1000<- length(intersect(topTab$ID, rownames(ranked)[1:1000]))
c100 <- length(intersect(topTab$ID[1:100], rownames(ranked)[1:100]))
c25 <- length(intersect(topTab$ID[1:25], rownames(ranked)[1:25]))


## ----loadAllData---------------------------------------------------------------
stopifnot(require(Biobase)) #library(Biobase)
load(file=file.path(dataDir, "normalized.filtered.Rda"))
my.eset <- eset_filtered


## ----analyzeAll1---------------------------------------------------------------
stopifnot(require(limma))
#designMatrix
lev<-as.factor(paste(pData(my.eset)$treat, pData(my.eset)$age, sep="."))
design <-model.matrix(~0+lev)
colnames(design)<-as.character(unique(lev))
rownames(design) <-rownames(pData(my.eset))
print(design)
#contrastsMatrix
cont.matrix <- makeContrasts (LPSvsMED.Aged=LPS.Aged-MED.Aged,
                              LPSvsMED.Young=LPS.Young-MED.Young,
                              AgedvsYoung.MED = MED.Aged-MED.Young,
                              levels=design)
cont.matrix



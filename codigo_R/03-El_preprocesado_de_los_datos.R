## ----losDatos, results='hide'-----------------------------------------------------------------------------------------------------------
if(!require(BiocManager)) install.packages("BiocManager")
if (!(require(CCl4))){
 BiocManager::install("CCl4")
}
if (!(require(estrogen))){
 BiocManager::install("estrogen")
}
if (!(require(affy))){
 BiocManager::install("affy")
}
if (!(require(affyPLM))){
 BiocManager::install("affyPLM")
}


## ----leeDatos1colors, results='hide', echo=TRUE, message=FALSE--------------------------------------------------------------------------
library(estrogen)
library(affy)
affyPath <- system.file("extdata", package = "estrogen")
adfAffy = read.AnnotatedDataFrame("phenoData.txt", sep="",  path=affyPath)
affyTargets = pData(adfAffy)
affyTargets$filename = file.path(affyPath, row.names(affyTargets))
affyRaw <- read.affybatch(affyTargets$filename, phenoData=adfAffy)
# show(affyRaw)
actualPath <- getwd()
setwd(affyPath)
allAffyRaw <- ReadAffy()
setwd(actualPath)


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------------------
class(affyRaw)
print(affyRaw)


## ----leeDatos2colores, results='hide', echo=TRUE, message=FALSE-------------------------------------------------------------------------
library("limma")
library("CCl4")
dataPath = system.file("extdata", package="CCl4")
adf = read.AnnotatedDataFrame("samplesInfo.txt", 
    path=dataPath)
#adf
targets = pData(adf)
targets$FileName = row.names(targets)
RG <- read.maimages(targets, path=dataPath, source="genepix")
attach(RG$targets)
newNames <-paste(substr(Cy3,1,3),substr(Cy5,1,3),substr(FileName,10,12), sep="")
colnames(RG$R)<-colnames(RG$G)<-colnames(RG$Rb)<-colnames(RG$Gb)<-rownames(RG$targets)<- newNames
# show(RG)


## ---------------------------------------------------------------------------------------------------------------------------------------
class(RG)
print(RG)


## ----plotAffyHist,echo=TRUE-------------------------------------------------------------------------------------------------------------
affySampleNames <- rownames(pData(allAffyRaw))
affyColores <- c(1,2,2,3,3,4,4,8,8)
affyLineas <- c(1,2,2,2,2,3,3,3,3)
hist(allAffyRaw, main="Signal distribution", col=affyColores, lty=affyLineas)
legend (x="topright", legend=affySampleNames , col=affyColores, lty=affyLineas, cex=0.7)


## ----plotAffyBoxplot--------------------------------------------------------------------------------------------------------------------
boxplot(allAffyRaw, main="Signal distribution", col=affyColores, las=2)


## ----plotPCAdef, results='hide'---------------------------------------------------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),
       ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10),
       )
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=0.8)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}


## ----tempNorm, echo= TRUE, results='hide'-----------------------------------------------------------------------------------------------
if (!file.exists("datos/affyNorm")){
  allAffyNorm<- rma(allAffyRaw)
  affyNorm <- rma(affyRaw)
  save(allAffyNorm, affyNorm, file="datos/affyNorm.Rda")
}else{
  load(file="datos/affyNorm.Rda")
}


## ----plotPCA2D1, echo=TRUE--------------------------------------------------------------------------------------------------------------
plotPCA(exprs(allAffyNorm), labels=affySampleNames, dataDesc="Grafico con todos los arrays\nIncluye muestra defectuosa")


## ----plotPCA2D2, echo=TRUE--------------------------------------------------------------------------------------------------------------
plotPCA(exprs(affyNorm), labels=colnames(exprs(affyNorm)), dataDesc="Grafico con todos los arrays\nSin la muestra defectuosa")


## ----plotDendro, echo=TRUE, fig.cap="Un cluster jerárquico sirve para determinar si las muestras se agrupan de forma natural según los grupos experimentales o si lo hacen de otra forma"----
clust.euclid.average <- hclust(dist(t(exprs(affyNorm))),method="average")
plot(clust.euclid.average, labels=colnames(exprs(affyNorm)), main="Hierarchical clustering of samples",  hang=-1, cex=0.7)


## ----c05normScatterMA, fig.cap="(a) Gráfico de  un canal frente al otro  (b) MA-plot (intensidad frente log-ratio)"---------------------
R <- RG$R[,"DMSCCl319"]
G <- RG$G[,"DMSCCl319"]
logR <- log(R)
logG <- log(G)
M <- logR-logG
A <- 0.5*(logR+logG)
opt<- par(mfrow=c(1,2))
plot(logR~ logG, main= "Scatterplot")
abline(lm(logR~logG), col="yellow")
plot(M~A, main= "MA-Plot")
abline(h=0, col="yellow")
par(opt)


## ----c05normScatterMA2, fig.cap="(a) Gráfico de  un canal frente al otro  (b) MA-plot (intensidad frente log-ratio). Estos gráficos sugieren la necesidad de normalizar los datos"----
library(marray)
data(swirl)
R <- swirl@maRf[,1]
G <- swirl@maGf[,1]
logR <- log(R)
logG <- log(G)
M <- logR-logG
A <- 0.5*(logR+logG)
opt<- par(mfrow=c(1,2))
plot(logR~ logG, main= "Scatterplot")
abline(lm(logR~logG), col="yellow")
plot(M~A, main= "MA-Plot")
abline(h=0, col="yellow")
par(opt)


## ----c05signalNoise, fig.cap="La relación señal ruido sirve para detectar posibles anormalidades o un background excesivamente alto como medida de calidad", echo=FALSE----
knitr::include_graphics("figures/c05signalNoise.png")


## ----c05plotAffy, fig.cap="Imágenes de cuatro microarrays de Affymetrix", echo=FALSE----------------------------------------------------
knitr::include_graphics("figures/c05plotAffy.png")


## ----c05MAPlotAffy, fig.cap="En los chips de Affymetrix la única forma de definir M (el log ratio) es comparar entre diferentes arrays", echo=FALSE----
knitr::include_graphics("figures/c05MAPlotAffy.png")


## ----affyPLM, echo=TRUE, results='hide'-------------------------------------------------------------------------------------------------
stopifnot(require(affyPLM))
Pset<- fitPLM(allAffyRaw)


## ----plotPLM ,fig.cap ="Graficos de diagnóstico calculados a nivel de sondas PLM", echo=TRUE--------------------------------------------
opt<- par(mfrow=c(2,1))
RLE(Pset, main = "Relative Log Expression (RLE)", 
    names=rownames(pData(allAffyRaw)), las=2, cex.axis=0.6)
NUSE(Pset, main = "Normalized Unscaled Standard Errors (NUSE)", las=2, 
     names=rownames(pData(allAffyRaw)), las=2, cex.axis=0.6)
par(opt)


## ----c06quantilNorm, fig.cap="El método RMA incluye una normalización por cuantiles como la representada esquemáticamente en esta figura", echo=FALSE----
knitr::include_graphics("figures/c06quantilNorm.png")


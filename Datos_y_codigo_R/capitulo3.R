### R code from vignette source 'capitulo6.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: losDatos (ACTUALIZADO)
###################################################
if (!require(BiocManager))  install.packages("BiocManager")
if (!require(CCl4)) BiocManager::install("CCl4")
if (!require(estrogen)) BiocManager::install("estrogen")
if (!require(affy)) BiocManager::install("affy")
if (!require(marray)) BiocManager::install("marray")
if (!require(limma)) BiocManager::install("limma")



###################################################
### code chunk number 2: leeDatos1colors
###################################################

require(estrogen)
require(affy)
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


###################################################
### code chunk number 3: leeDatos2colores
###################################################
require("limma")
require("CCl4")
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


###################################################
### code chunk number 4: plotHist
###################################################
affySampleNames <- rownames(pData(allAffyRaw))
affyColores <- c(1,2,2,3,3,4,4,8,8)
affyLineas <- c(1,2,2,2,2,3,3,3,3)
hist(allAffyRaw, main="Signal distribution", col=affyColores, lty=affyLineas)
legend (x="topright", legend=affySampleNames , col=affyColores, lty=affyLineas, cex=0.7)


###################################################
### code chunk number 5: plotBoxplot
###################################################
boxplot(allAffyRaw, main="Signal distribution", col=affyColores, las=2)


###################################################
### code chunk number 6: plotPCAdef
###################################################
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


###################################################
### code chunk number 7: tempNorm
###################################################
if (!file.exists("allAffyNorm")){
  allAffyNorm<- rma(allAffyRaw)
  affyNorm <- rma(affyRaw)
  save(allAffyNorm, affyNorm, file="affyNorm.Rda")
}else{
  load(file="affyNorm.Rda")
}


###################################################
### code chunk number 8: plotPCA2D
###################################################
opt <- par(mfrow=c(2,1))
plotPCA(exprs(allAffyNorm), labels=affySampleNames, dataDesc="PCA for all arrays\nincludes defective sample")
plotPCA(exprs(affyNorm), labels=colnames(exprs(affyNorm)), dataDesc="PCA for all arrays")
par(opt)


###################################################
### code chunk number 9: plotDendro
###################################################
clust.euclid.average <- hclust(dist(t(exprs(affyNorm))),method="average")
plot(clust.euclid.average, labels=colnames(exprs(affyNorm)), main="Hierarchical clustering of samples",  hang=-1, cex=0.7)


###################################################
### code chunk number 10: plotDeg (ACTUALIZADO)
###################################################

### Esta opcion se ha desactivado, puesto que parece dar problemas.
### El grafico de degradación ha caído en desuso porque los arrays modernos tienen sondas en todos o en la mayoría de los exones
#
# deg<-AffyRNAdeg(allAffyRaw, log.it=T)
# colores<-c("red", rep("blue",8))
# lineas <- c(1, rep(3,8))
# plotAffyRNAdeg(deg, col=colores, lty=lineas)
# legend (x="bottomleft", legend=rownames(pData(allAffyRaw)), col=colores, lty=lineas, cex=0.7)
# 

###################################################
### code chunk number 11: affyPLM
###################################################
stopifnot(require(affyPLM))
Pset<- fitPLM(allAffyRaw)


###################################################
### code chunk number 12: plotPLM (ACTUALIZADO)
###################################################

opt<- par(mfrow=c(2,1))
RLE(Pset, main = "Relative Log Expression (RLE)", 
    names=rownames(pData(allAffyRaw)), las=2, cex.axis=0.6)
NUSE(Pset, main = "Normalized Unscaled Standard Errors (NUSE)",  
     names=rownames(pData(allAffyRaw)), las=2, cex.axis=0.6)
par(opt)

###############################################################
### code chunk number 1: marrayPlots.Rnw:105-108 (ACTUALIZADO)
###############################################################

library(marray)
data(swirl) 
maPlate(swirl)<-maCompPlate(swirl,n=384) 


###################################################
### code chunk number 2: marrayPlots.Rnw:128-131
###################################################

Gcol<- maPalette(low="white", high="green",k=50) 
Rcol<- maPalette(low="white", high="red", k=50) 
RGcol<-maPalette(low="green", high="red", k=50) 


###################################################
### code chunk number 3: maImageGb
###################################################

tmp<-image(swirl[,3], xvar="maGb", subset=TRUE, col=Gcol,contours=FALSE, bar=FALSE) 


###################################################
### code chunk number 4: maImageRb
###################################################
tmp<-image(swirl[,3], xvar="maRb", subset=TRUE, col=Rcol, contours=FALSE, bar=FALSE) 


###################################################
### code chunk number 5: maImageMraw1
###################################################
tmp<-image(swirl[,3], xvar="maM", bar=FALSE, main="Swirl array 93: image of pre--normalization M") 


###################################################
### code chunk number 6: maImageMraw2
###################################################
tmp<-image(swirl[,3], xvar="maM", subset=maTop(maM(swirl[,3]), h=0.10,
l=0.10), col=RGcol, contours=FALSE, bar=FALSE,main="Swirl array 93:
image of pre--normalization M for  10 %  tails")  


###################################################
### code chunk number 7: maImageSpotCol
###################################################
tmp<- image(swirl[,3], xvar="maSpotCol", bar=FALSE) 


###################################################
### code chunk number 8: maImagePrintTip
###################################################
tmp<- image(swirl[,3], xvar="maPrintTip", bar=FALSE) 


###################################################
### code chunk number 9: maImageControls
###################################################
tmp<- image(swirl[,3], xvar="maControls",col=heat.colors(10),bar=FALSE) 


###################################################
### code chunk number 10: maImagePlate
###################################################
tmp<- image(swirl[,3], xvar="maPlate",bar=FALSE) 


###################################################
### code chunk number 11: maBoxplot1pre
###################################################
boxplot(swirl[,3], xvar="maPrintTip", yvar="maM", main="Swirl array 93: pre--normalization") 


###################################################
### code chunk number 12: maBoxplot2pre
###################################################
boxplot(swirl, yvar="maM", main="Swirl arrays: pre--normalization") 


###################################################
### code chunk number 13: marrayPlots.Rnw:343-344
###################################################
swirl.norm <- maNorm(swirl, norm="p")


###################################################
### code chunk number 14: maBoxplot1post
###################################################
boxplot(swirl.norm[,3], xvar="maPrintTip", yvar="maM",
	main="Swirl array 93: post--normalization") 

opt<- par(mfrow=c(1,2))
boxplot(swirl, yvar="maM", main="Swirl arrays: pre--normalization") 
boxplot(swirl.norm, yvar="maM", col="green", main="Swirl arrays: post--normalization")
par(opt)

###################################################
### code chunk number 15: maBoxplot2post
###################################################
par(mfrow=c(1,1))
boxplot(swirl.norm, yvar="maM", col="green", main="Swirl arrays: post--normalization") 


###################################################
### code chunk number 16: maPlot1pre
###################################################
defs<-maDefaultPar(swirl[,3],x="maA",y="maM",z="maPrintTip")

# Function for plotting the legend
legend.func<-do.call("maLegendLines",defs$def.legend)

# Function for performing and plotting lowess fits
lines.func<-do.call("maLowessLines",c(list(TRUE,f=0.3),defs$def.lines))

plot(swirl[,3], xvar="maA", yvar="maM", zvar="maPrintTip",
		      lines.func,
		      text.func=maText(),
		      legend.func,
		      main="Swirl array 93: pre--normalization MA--plot") 


###################################################
### code chunk number 17: maPlot1post
###################################################
plot(swirl.norm[,3], xvar="maA", yvar="maM", zvar="maPrintTip",
		      lines.func,
		      text.func=maText(),
		      legend.func,
		      main="Swirl array 93: post--normalization MA--plot") 



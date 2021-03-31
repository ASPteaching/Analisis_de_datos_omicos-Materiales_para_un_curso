## ----loadValues-----------------------------------------------------
require(Biobase)
require(limma)
load("datos/celltypes/celltypes-normalized.rma.Rda")
load("datos/celltypes/celltypes-fit.main.Rda")
topTab_LPS.in.AGED <- topTable (fit.main, number=nrow(fit.main), coef="LPS.in.AGED", adjust="fdr",lfc=2)
topTab_LPS.in.YOUNG <- topTable (fit.main, number=nrow(fit.main), coef="LPS.in.YOUNG", adjust="fdr",lfc=2)
topTab_AGE  <- topTable (fit.main, number=nrow(fit.main) , coef="AGE", adjust="fdr", lfc=2)


## ----decideTests01, echo=F------------------------------------------
require(limma)
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.0001, lfc=1)


## ----resumeDecideTests1---------------------------------------------
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,]
print(summary(res))


## ----decideTests02, echo=F, eval=TRUE-------------------------------
res<-decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.0001, lfc=2)
sum.res.rows<-apply(abs(res),1,sum)
res.selected<-res[sum.res.rows!=0,]
print(summary(res))


## ----venn01---------------------------------------------------------
vennDiagram (res.selected[,1:3], main="Genes in common #1", cex=0.9)


## ----anota01, print=FALSE-------------------------------------------
library(mouse4302.db)
anotData <- capture.output(mouse4302())
print(anotData)
cat ("... output continues until ", length(anotData), " lines.\n")


## ----anota02, print=FALSE-------------------------------------------
if (!(require(org.Mm.eg.db))){
        biocLite("org.Mm.eg.db")
      }
require(org.Mm.eg.db)
anotData <- capture.output(org.Mm.eg())
print(anotData[1:15])
cat ("... output continues until ", length(anotData), " lines.\n")


## ----top5Genes01----------------------------------------------------
top5 <-topTab_LPS.in.AGED$ID[1:5]
cat("Usando mget\n")
geneSymbol5.1 <- unlist(mget(top5, mouse4302SYMBOL))
geneSymbol5.1
cat("Usando toTable\n")
genesTable<- toTable(mouse4302SYMBOL)
rownames(genesTable) <-  genesTable$probe_id
genesTable[top5, 2]
cat("Usando getSYMBOL\n")
require(annotate)
geneSymbol5.3 <- getSYMBOL(top5, "mouse4302.db")
geneSymbol5.3


## ----annaffy01, eval=FALSE------------------------------------------
## require(annaffy)
## genesSelected <- rownames(res.selected)
## at <- aafTableAnn(genesSelected, "mouse4302.db")
## saveHTML (at, file="results/anotations.html",
##           "Annotations for selected genes")


## ----prepareData01--------------------------------------------------
probeNames<-rownames(res)
probeNames.selected<-probeNames[sum.res.rows!=0]
exprs2cluster <-exprs(eset_rma_filtered)[probeNames.selected,]
colnames(exprs2cluster)<- c("OldLPS80L", "OldLPS86L", "OldLPS88L",
                            "OldMED81m", "OldMED82m", "OldMED84m",
                            "YouLPS75L",  "YouLPS76L",  "YouLPS77L",
                            "YouMED71m", "YouMED72m", "YouMED73m")


## ----plotHeatMap01, fig=T-------------------------------------------
color.map <- function(grupo) {
  switch(grupo,
         "yellow",
         "red",
         "blue",
         "pink")
}
grupColors <- unlist(lapply(pData(eset_rma_filtered)$grupo, color.map))
heatmap(exprs2cluster,
        cexCol=0.8,
        main="Heatmap para las tres comparaciones de 'celltypes'", cex.main=0.8)


## ----plotHeatMap02, fig=T-------------------------------------------
require("gplots")
heatmap.2(exprs2cluster,
          col=bluered(75), scale="row",
          ColSideColors=grupColors, key=TRUE, symkey=FALSE,
          density.info="none", trace="none", cexCol=0.8,
          main="Heatmap de las muestras de 'celltypes'",   cex.main=0.6)


## ----GOanots01------------------------------------------------------
require(annotate)
(top1 <-rownames(topTab_LPS.in.AGED)[1])
(geneSymbol1 <- getSYMBOL(top1, "mouse4302.db"))
GOAnots1 <- mget(top1, mouse4302GO)
for (i in 1:length(GOAnots1)){
  for (j in 1:length(GOAnots1[[i]])){
    GOAnot <- GOAnots1[[i]][[j]][[1]]
    cat(top1[i],geneSymbol1[i],GOAnot,substr(Term(GOAnot),1,30), "\n")
  }
}


## ----GOAnalysis01, eval=FALSE---------------------------------------
## require(GOstats)
## require(mouse4302.db)
## require(org.Mm.eg.db)
## 
##   # Seleccionamos la "topTable"
##   topTab <- topTab_LPS.in.AGED
##   # Definimos el universo de genes: todos los que se han incluido en el análisis
##   # EL programa trabaja con identificadores "entrez" y no admite duplicados
## 
##   entrezUniverse <- unique(getEG(as.character(rownames(topTab)), "mouse4302.db"))
## 
##   # Escogemos los grupos de sondas a incluir en el análisis
##   # Este análisis trabaja bien con varios centenares de genes
##   # por lo que es habitual basarse en p-valores sin ajustar para incluirlos
## 
##   whichGenes<-topTab["adj.P.Val"]<0.001
##   geneIds <-   unique(getEG(as.character(rownames(topTab)[whichGenes]),"mouse4302.db"))
## 
##   # Creamos los "hiperparámetros" en que se basa el análisis
##   GOparams = new("GOHyperGParams",
##     geneIds=geneIds, universeGeneIds=entrezUniverse,
##     annotation="org.Mm.eg.db", ontology="BP",
##     pvalueCutoff=0.001, conditional=FALSE,
##     testDirection="over")
## 
##   # Ejecutamos los análisis
## 
##   GOhyper = hyperGTest(GOparams)
## 
## # Creamos un informe html con los resultados
##    comparison = "topTab_LPS.in.AGED"
##    GOfilename =file.path(resultsDir,
##      paste("GOResults.",comparison,".html", sep=""))
##   htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))


## ----KEGGAnalysis01, eval=FALSE-------------------------------------
## 
##   KEGGparams = new("KEGGHyperGParams",
##     geneIds=geneIds, universeGeneIds=entrezUniverse,
##     annotation="org.Mm.eg.db",
##     pvalueCutoff=0.01, testDirection="over")
## 
##   # Ejecutamos los análisis
## 
##   KEGGhyper = hyperGTest(KEGGparams)
## 
## # Creamos un informe html con los resultados
##  comparison = "topTab_LPS.in.AGED"
##  KEGGfilename =file.path(resultsDir,
##      paste("KEGGResults.",comparison,".html", sep=""))
##   htmlReport(KEGGhyper, file = KEGGfilename, summary.args=list("htmlLinks"=TRUE))


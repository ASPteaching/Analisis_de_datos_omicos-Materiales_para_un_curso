
installifnot <- function (packageName){
  if (!(require(packageName, character.only=TRUE))) {
    install.packages(packageName, dep=TRUE)
  }else{
    detach(paste ("package", packageName, sep=":"), character.only=TRUE)
  } 
}
bioCifnot <- function (packageName){
  if (!(require(packageName, character.only=TRUE))) {
  BiocManager::install(packageName)
  }else{
    detach(paste ("package", packageName, sep=":"), character.only=TRUE)
  }  
}
installifnot("knitr")
installifnot("ggrepel")
bioCifnot ("affy")
bioCifnot ("genefilter")
bioCifnot ("multtest")
bioCifnot ("arrayQualityMetrics")
bioCifnot ("limma")
bioCifnot ("estrogen")
bioCifnot ("CCl4")
bioCifnot ("annotate")
bioCifnot ("annaffy")
bioCifnot ("org.Hs.eg.db")
bioCifnot ("org.Mm.eg.db")
bioCifnot ("hgu95av2.db")
bioCifnot ("hgu133a.db")
bioCifnot ("mouse4302.db")
bioCifnot ("GO.db")
bioCifnot ("KEGG.db")
bioCifnot ("Reactome.db")
bioCifnot ("GOstats")
```
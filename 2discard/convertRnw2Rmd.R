curDir<- getwd()
library(Hmisc)
getRs('convertL2M.r', put='source')   # gets from github
# setwd("/media/alex/Seagate Expansion Drive/Treballs/2011-Materials en elaboracio per assignatura ANALISI DE DADES D'ALT RENDIMENT/1.1. Introduccio Als microarrays-c1")
convertL2M('04-Seleccion_de_genes_DE.Rnw', '04-Seleccion_de_genes_DE.Rmd')
setwd(curDir)

library(stringr)

fesCanvis <- function (inputFile, outputFile){
  x <- readLines(inputFile)
  y1 <- str_replace_all(x, "\\\\'a", "á")
  y2 <- str_replace_all(y1, "\\\\'e", "é")
  y3 <- str_replace_all(y2, "\\\\'i", "í")
  y4 <- str_replace_all(y3, "\\\\'o", "ó")
  y5 <- str_replace_all(y4, "\\\\'u", "ú")
  # y6 <- str_replace_all(y5, "\\\~n", "ñ")
  # y6 <- str_replace_all(y6, "\\\\`\\\\`", '"')
  
  cat(y4, file=outputFile, sep="\n")
}

fesCanvis("capitulo7.rnw", "04-Seleccion_de_genes_DE.Rmd")

##########################33

x1 <- "cap\'itulo"
x2 <- "\\section{Fuentes de variabilidad}"
pattern <- "\\section\\{[[:alnum:][:blank:]]+\\}"
replacement <-"\\s[[:alnum:][:blank:]]+"
str_replace_all( x2, pattern, replacement)


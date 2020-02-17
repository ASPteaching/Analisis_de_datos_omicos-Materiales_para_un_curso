# Omics_Data_Analysis-Course_Materials

## Presentation

This repository contains R code and materials for my -never published- book on Microarray Data Analysis.

It will be progressively filled with 

- R scripts for those chapters where there is R code
- Rmarkdown files to generate the book

The book is almost ten years old, what, in this field, means it is, at least partly, outdated. 

- The first thing I will do is to update the scripts so that they work with 2019's Bioconductor. 
- After this I will re-write the old .rnw (Sweave) files and create a bookdown document.

## Errata detected and changes applied

### ID vs rownames
- Chapter 5 and 6 contained references to a field named "ID" in a topTable (indeed in the "lmFit" object from where the topTable was extracted.)
- This was so in old limma versions. A few years ago limma maintainers decided to remove this field and keep the probeset identifiers as rownames.
- I have updated the code to remove references to this "ID" and change them by "rownames()"


require(Hmisc)
getRs('convertL2M.r', put='source')   # gets from github
convertL2M('Practical-Bioconductor-Classes.Rnw', 'Practical-Bioconductor-Classes.Rmd')

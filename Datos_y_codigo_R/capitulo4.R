### R code from vignette source 'capitulo7.rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: c7codec01 (ACTUALIZADO con un "supressWarnings")
###################################################

stopifnot(require(Biobase)) 
load ("datos/celltypes-normalized.rma.Rda")
suppressWarnings(my.eset <- eset_rma_filtered)
grupo_1 <- as.factor(pData(my.eset)$treat)
stopifnot(require(genefilter))
teststat <-rowttests(my.eset, "treat")
print(teststat[1:5,])


###################################################
### code chunk number 2: sortBypvals
###################################################

ranked <-teststat[order(teststat$p.value),]
print(ranked[1:5,])


###################################################
### code chunk number 3: selectedGenes
###################################################

selectedTeststat <- ranked[ranked$p.value < 0.01,]


###################################################
### code chunk number 4: adjustPvals
###################################################
stopifnot(require(multtest))
procs <- c("Bonferroni","BH", "BY")
adjPvalues <- mt.rawp2adjp(teststat$p.value, procs)
names(adjPvalues)
ranked.adjusted<-cbind(ranked, adjPvalues$adjp)
head(ranked.adjusted)


###################################################
### code chunk number 5: selectAdjusted
###################################################

selectedAdjusted<-ranked.adjusted[ranked.adjusted$BY<0.001,]
nrow(selectedAdjusted)


###################################################
### code chunk number 6: selectLimma
###################################################

require(Biobase)
require(limma)
load ("datos/celltypes-normalized.rma.Rda")
my.eset <- eset_rma_filtered
targets <- pData(my.eset)
age.treat<- paste(targets$age,targets$treat, sep=".")
lev<-factor(age.treat, levels=unique(age.treat))
design <-model.matrix(~0+lev)
colnames(design)<-levels(lev)
rownames(design) <-rownames(targets)
print(design)


###################################################
### code chunk number 7: setContrasts
###################################################

require(limma)
cont.matrix <- makeContrasts (
      LPS.in.AGED=(Aged.LPS-Aged.MED),
      LPS.in.YOUNG=(Young.LPS-Young.MED),
      AGE=(Aged.MED-Young.MED),
      levels=design)
cont.matrix


###################################################
### code chunk number 8: linearmodelfit
###################################################

require(limma)
fit<-lmFit(my.eset, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)
save(fit.main, file="results/celltypes-fit.main.Rda")


###################################################
### code chunk number 9: capitulo7.rnw:1132-1135
###################################################

topTab_LPS.in.AGED <- topTable (fit.main, number=nrow(fit.main), coef="LPS.in.AGED", adjust="fdr")
topTab_LPS.in.YOUNG <- topTable (fit.main, number=nrow(fit.main), coef="LPS.in.YOUNG", adjust="fdr")
topTab_AGE  <- topTable (fit.main, number=nrow(fit.main) , coef="AGE", adjust="fdr")

###################################################
### code chunk number 10: volcano1
###################################################

coefnum = 1
opt <- par(cex.lab = 0.7)
volcanoplot(fit.main, coef=coefnum, highlight=10, names=fit.main$ID, 
            main=paste("Differentially expressed genes",colnames(cont.matrix)[coefnum], sep="\n"))
abline(v=c(-1,1))
par(opt)


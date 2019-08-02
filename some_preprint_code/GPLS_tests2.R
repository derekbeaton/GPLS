### this all sort of works...
## just ever so different in slight ways from expected.

rm(list=ls())
library(GSVD)
library(ExPosition)
library(ours)
source('plsrca.fin.r')
source('gplss.R')

load("../../ADNI_DATA/TADPOLE.fin2.rda")
load("../../ADNI_DATA/genetic.data.rda")
load("../../ADNI_DATA/TADPOLE_Dict.fin2.rda")


genetic.data.orig <- genetic.data
genetic.data <- genetic.data[,-1]

ma <- tapply(genetic.data,col(genetic.data),function(x){table(unlist(strsplit(as.vector(x),"")))},simplify=F)	
names(ma) <- colnames(genetic.data)		
ma.finder <- lapply(ma,function(x){ a<-x/sum(x); a==min(a) })
geno.map <- do.call(rbind,lapply(ma.finder,
                                 function(x){ 
                                   c( 
                                     paste(rep(names(x)[which(x==T)],2),collapse=""),
                                     paste(c(names(x)[which(x==T)],names(x)[which(x==F)]),collapse=""),
                                     paste(rep(names(x)[which(x==F)],2),collapse="")
                                   )
                                 }))
colnames(geno.map) <- c("aa","Aa","AA")

maf <- .05
find.low.genos <- tapply(genetic.data,col(genetic.data),function(x){ summary(as.factor(x[!is.na(x)])) < (length(x)*maf) },simplify=F)
names(find.low.genos) <- colnames(genetic.data)

low.geno.locs <- which(unlist(lapply(find.low.genos,function(x){any(x)})))

for(i in 1:length(low.geno.locs)){
  
  genetic.data[,low.geno.locs[i]] <- ifelse(
    genetic.data[,low.geno.locs[i]] == geno.map[names(low.geno.locs[i]),"AA"],
    geno.map[names(low.geno.locs[i]),"AA"],
    paste(geno.map[names(low.geno.locs[i]),"Aa"],geno.map[names(low.geno.locs[i]),"aa"],sep="+")
  )
  
}

genetic.nom <- make.data.nominal(genetic.data)
genetic.nom.design <- make.data.nominal(as.matrix(attributes(genetic.nom)$variable.map))
rownames(genetic.nom.design) <- colnames(genetic.nom)

VENT.FRAME <- data.frame(VENT=TADPOLE.fin[,101],ICV=TADPOLE.fin[,104])
HIPPO.FRAME <- data.frame(HIPPO=TADPOLE.fin[,102],ICV=TADPOLE.fin[,104])
WB.FRAME <- data.frame(WB=TADPOLE.fin[,103],ICV=TADPOLE.fin[,104])
rownames(VENT.FRAME) <- rownames(HIPPO.FRAME) <- rownames(WB.FRAME) <- rownames(TADPOLE.fin)
VENT <- resid(lm(VENT~ICV,VENT.FRAME)) + mean(TADPOLE.fin[,101],na.rm=T)
HIPPO <- resid(lm(HIPPO~ICV,HIPPO.FRAME)) + mean(TADPOLE.fin[,102],na.rm=T)
WB <- resid(lm(WB~ICV,WB.FRAME)) + mean(TADPOLE.fin[,103],na.rm=T)
structural.norm <- matrix(NA,nrow(TADPOLE.fin),3)
colnames(structural.norm) <- c("VENT","HIPPO","WB")
rownames(structural.norm) <- rownames(TADPOLE.fin)
structural.norm[names(VENT),"VENT"] <- VENT
structural.norm[names(HIPPO),"HIPPO"] <- HIPPO
structural.norm[names(WB),"WB"] <- WB

## HARDCODE HACK.
grp.cols <- (viridisLite::cividis(20))[seq(1,20,4)+2][c(3,5,1,4,2)]
names(grp.cols) <- unique(TADPOLE.fin$DX_bl)
ind.cols <- plyr::mapvalues(TADPOLE.fin$DX_bl,from=names(grp.cols),to=grp.cols)
names(ind.cols) <- rownames(TADPOLE.fin)

DESIGN <- make.data.nominal(as.matrix(TADPOLE.fin$DX_bl))
colnames(DESIGN) <- gsub("\\.","",colnames(DESIGN))
W.DESIGN <- apply(DESIGN,2,function(x){x/sum(x)})
rownames(W.DESIGN) <- rownames(DESIGN) <- rownames(TADPOLE.fin)



CTX.uptake <- TADPOLE.fin[,grep("CTX",colnames(TADPOLE.fin))]
colnames(CTX.uptake) <- gsub("_UCBERKELEYAV45_10_17_16","",gsub("CTX_","",colnames(CTX.uptake)))

drop.rows <- which(rowSums(is.na(CTX.uptake))==70)
drop.cols <- grep("UNKNOWN",colnames(CTX.uptake))
CTX.uptake_drop <- CTX.uptake[-c(drop.rows),-c(drop.cols)]
CTX.uptake_drop <- apply(CTX.uptake_drop,2,function(x){ x[which(is.na(x))] <- mean(x,na.rm=T); x })

genetic.nom_drop <- genetic.nom[rownames(CTX.uptake_drop),]



age.doubled <- escofier.coding(as.matrix(TADPOLE.fin$AGE[-c(drop.rows)]),scale=T)
colnames(age.doubled) <- c('AGE-','AGE+')
edu.thermometer <- thermometer.coding(as.matrix(TADPOLE.fin$PTEDUCAT[-c(drop.rows)]),norm.to.one = T)
colnames(edu.thermometer) <- c('EDU+','EDU-')
beh.brain.confounds <- cbind(makeNominalData(as.matrix(TADPOLE.fin[-c(drop.rows),c("PTGENDER")])),age.doubled,edu.thermometer)
colnames(beh.brain.confounds) <- gsub("\\.","",colnames(beh.brain.confounds))
rownames(beh.brain.confounds) <- rownames(TADPOLE.fin)[-c(drop.rows)]




Y <- CTX.uptake_drop
X <- beh.brain.confounds

ZZ.res <- ca.preproc(genetic.nom_drop)


plsrca.res <- plsrca(X,Y)

### the clue is here... Zx & Zy are not equivalent here to the standard Zs.
X.res <- ca.preproc(X)
Y.res <- ca.preproc(Y)
ZX.in <- sweep(X.res$Zx, 1, sqrt(X.res$m), "/")
ZY.in <- sweep(Y.res$Zx, 1, sqrt(Y.res$m), "/")

Zx.in.test <- (diag(X.res$m) %^% (-1/2)) %*% X.res$Zx
Zy.in.test <- (diag(Y.res$m) %^% (-1/2)) %*% Y.res$Zx

gplsrca.res <- gplsreg(ZX.in, ZY.in, 1/X.res$w, 1/Y.res$w, comps = length(plsrca.res$d.vec), tol = .Machine$double.eps)
# gplsrca.res_can <- gplscan(ZX.in, ZY.in, 1/X.res$w, 1/Y.res$w, comps = length(plsrca.res$d.vec), tol = .Machine$double.eps)
# 
# crossprod(gplsrca.res_can$LX)
# crossprod(gplsrca.res_can$LY)
# diag(t(gplsrca.res_can$LX) %*% gplsrca.res_can$LY )
# gplsrca.res_can$d


gplsrca.res$fi / plsrca.res$fi.matrix
gplsrca.res$fj / plsrca.res$fj.matrix
gplsrca.res$u / plsrca.res$u.matrix
gplsrca.res$p / plsrca.res$p.matrix
gplsrca.res$v / plsrca.res$v.matrix
gplsrca.res$q / plsrca.res$q.matrix
gplsrca.res$LX / plsrca.res$lx
gplsrca.res$LY / plsrca.res$ly
gplsrca.res$d / plsrca.res$d.vec
gplsrca.res$b / plsrca.res$beta.vec
gplsrca.res$pred.U / plsrca.res$pred.u.matrix
gplsrca.res$r2.x.cumulative / plsrca.res$r2.x.cumulative
gplsrca.res$r2.y.cumulative / plsrca.res$r2.y.cumulative
gplsrca.res$Tmat / plsrca.res$Xt.matrix
gplsrca.res$Y.recs / plsrca.res$Y.recs
gplsrca.res$Y.resids / plsrca.res$Y.resids
gplsrca.res$X.recs / plsrca.res$X.recs
gplsrca.res$X.resids / plsrca.res$X.resids
gplsrca.res$Y.rec / plsrca.res$Y.rec

  ## some reconstructions that have to happen within the plscar code not gplsreg
((gplsrca.res$Y.rec + Y.res$Ex) * sum(Y)) / plsrca.res$Y.hat
  ## no...
#gplsrca.res$Y.resid / plsrca.res$Y.resid
  ## yes
(((gplsrca.res$Y.resid/sqrt(nrow(Y)) + Y.res$Ex)) * sum(Y)) / plsrca.res$Y.resid


### quick test of gplsrca...

Yhat <- ((gplsrca.res$Y.rec + Y.res$Ex) * sum(Y))
Yresid <- (((gplsrca.res$Y.resid/sqrt(nrow(Y)) + Y.res$Ex)) * sum(Y))
Yresid2 <- ((gplsrca.res$Y.resid %*% diag(sqrt(1/Y.res$w)))  + Y.res$Ex) * sum(Y)




hat.ca <- ca(Yhat)
resid.ca <- ca(Yresid)
resid2.ca <- ca(Yresid2)
  ## near orthogonal
cor(hat.ca$fi,resid.ca$fi)


hat.ca2 <- ca(plsrca.res$Y.hat)
resid.ca2 <- ca(plsrca.res$Y.resid)
cor(hat.ca2$fi,resid.ca2$fi)


Xhat.ca2 <- ca(plsrca.res$X.hat)
Xresid.ca2 <- ca(plsrca.res$X.resid)
cor(Xhat.ca2$fi,Xresid.ca2$fi)


  ### NOTE: I need to force a stop of the 
  ## I don't know what this is...
rec.svd <- gsvd(gplsrca.res$Y.rec,RW=Y.res$w)
resid.svd <- gsvd(gplsrca.res$Y.resid,RW=Y.res$w)
  ## effectively orthogonal
cor(rec.svd$fi,resid.svd$fi)
  ### this needs to be noted in the "near orthogonal" part I mention...


### ok so rebuild gpls rec and resid based on Escofier...
new.Y.rec <- (gplsrca.res$Y.rec + Y.res$Ex) * sum(Y)
new.Y.resid <- (gplsrca.res$Y.resid + Y.res$Ex) * sum(Y)
  ### I should test this with X, too...

sum(new.Y.rec)/sum(Y)
rowSums(new.Y.rec)/rowSums(Y)
colSums(new.Y.rec)/colSums(Y)


sum(new.Y.resid)/sum(Y)
rowSums(new.Y.resid)/rowSums(Y)
colSums(new.Y.resid)/colSums(Y)


new.rec.ca <- ca(new.Y.rec)
new.resid.ca <- ca(new.Y.resid)

cor(new.rec.ca$u, new.resid.ca$u)
cor(new.rec.ca$fi, new.resid.ca$fi)

  ### should Y.rec + Y.resid = Y?
    ## what are the properties of Y.rec & Y.resid?
svd.rec <- gsvd(gplsrca.res$Y.rec)
svd.resid <- gsvd(gplsrca.res$Y.resid)
cor(svd.rec$u, svd.resid$u)
cor(svd.rec$fi, svd.resid$fi)




### plain CA with build back & deflation
ca.y <- ca(Y)

(ca.y$u %*% diag(ca.y$d) %*% t(ca.y$v))
(ca.y$u %*% diag(ca.y$d) %*% t(ca.y$v)) / Y.res$weightedZx
diag(sqrt(Y.res$m)) %*% (ca.y$u %*% diag(ca.y$d) %*% t(ca.y$v)) %*% diag(sqrt(Y.res$w))
(diag(sqrt(Y.res$m)) %*% (ca.y$u %*% diag(ca.y$d) %*% t(ca.y$v)) %*% diag(sqrt(Y.res$w))) / Y.res$Zx

(((diag(sqrt(Y.res$m)) %*% (ca.y$u %*% diag(ca.y$d) %*% t(ca.y$v)) %*% diag(sqrt(Y.res$w)))) + Y.res$Ex) / Y.res$Ox
((((diag(sqrt(Y.res$m)) %*% (ca.y$u %*% diag(ca.y$d) %*% t(ca.y$v)) %*% diag(sqrt(Y.res$w)))) + Y.res$Ex) * sum(Y)) / Y





(ca.y$p %*% diag(ca.y$d) %*% t(ca.y$q)) / Y.res$Zx
(ca.y$p %*% diag(ca.y$d) %*% t(ca.y$q)) / Y.res$weightedZx

(diag(sqrt(Y.res$m)) %*% (ca.y$p %*% diag(ca.y$d) %*% t(ca.y$q)) %*% diag(sqrt(1/Y.res$w))) / Y.res$weightedZx



ca.y$u / (diag(sqrt(Y.res$m)) %*% ca.y$p)
ca.y$v / (diag(sqrt(1/Y.res$w)) %*% ca.y$q)





  ## ok these are the same.
test.y.resid <- (plsrca.res$Zy.orig - plsrca.res$Y.rec)
test.y.rec <- (plsrca.res$Zy.orig - test.y.resid)

cor(test.y.resid, plsrca.res$Y.rec)
cor(plsrca.res$Zy.orig, plsrca.res$Y.rec)
cor(plsrca.res$Zy.orig, test.y.resid)

test.resid <- sweep(sweep(sweep(test.y.resid, 1, sqrt(plsrca.res$wiy), "*"),1,sqrt(1/plsrca.res$wiy),"*"),2,sqrt(1/plsrca.res$wjy),"*")
test.rec <- sweep(sweep(sweep(test.y.rec, 1, sqrt(plsrca.res$wiy), "*"),1,sqrt(1/plsrca.res$wiy),"*"),2,sqrt(1/plsrca.res$wjy),"*")


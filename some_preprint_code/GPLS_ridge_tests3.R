### this all sort of works...
## just ever so different in slight ways from expected.

rm(list=ls())
library(GSVD)
library(ExPosition)
library(ours)

source('plsrca.fin.r')
source('gplss.R')
source('rplsca.R')

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




# Y <- CTX.uptake_drop
# X <- beh.brain.confounds

# Y <- CTX.uptake_drop[,grep("LH_",colnames(CTX.uptake_drop))]
# X <- CTX.uptake_drop[,grep("RH_",colnames(CTX.uptake_drop))]

X <- beh.brain.confounds
Y <- genetic.nom_drop


X.res <- ca.preproc(X)
Y.res <- ca.preproc(Y)


#### ok first note: something is blowing up here where it shouldn't; find it and fix it
#### second note: square away various SVD issues (e.g., as set right now for plsrca() and takane's rplsrcar)
#### third note: clean this stuff up...

plsrca.res <- plsrca(X,Y) ### this might be fixed through the newest GSVD update but it looks like GSVD is out of date.
plsrca.res_flip <- plsrca(Y,X) ### this might be fixed through the newest GSVD update but it looks like GSVD is out of date.
plsrca.res_flip_tol <- plsrca(Y,X,tol = sqrt(.Machine$double.eps * 2)) ### this might be fixed through the newest GSVD update but it looks like GSVD is out of date.
plsrca.res_flip_inf <- plsrca(Y,X,tol = -Inf) ### this might be fixed through the newest GSVD update but it looks like GSVD is out of date.
  ### MAKE SURE TO TEST THIS WITH BEH.BRAIN AS Y and GENE as X.

### OK something still blowing up on $u, which means there may not be anything coming back...

gplsrca.res <- gplsreg( (diag(X.res$m)%^%(-1/2)) %*% X.res$Zx, (diag(Y.res$m)%^%(-1/2)) %*% Y.res$Zx, 1/X.res$w, 1/Y.res$w )

lambdas <- c(0, 1, 2, 3, 4, 5, 10, 25, 50)#, 75, 100, 200, 250, 500, 1000, 10000)

rplsres_all <- rplsca(X,Y,pls.type="reg",lambdas=lambdas,ridge.type="fast",comps=0, tol=.Machine$double.eps)
rplsres_all_takane <- rplsca(X,Y,pls.type="reg",lambdas=lambdas,ridge.type="",comps=0)
  ## sometimes this can't be done because of some Lapack failure for the SVD. Kind of out of my control.


## key comparisons
plsrca.res$u.matrix / rplsres_all$lambda_0$u
plsrca.res$v.matrix / rplsres_all$lambda_0$v
plsrca.res$d.vec / rplsres_all$lambda_0$d
plsrca.res$beta.vec / rplsres_all$lambda_0$b
plsrca.res$Xt.matrix / rplsres_all$lambda_0$Tmat
plsrca.res$lx / rplsres_all$lambda_0$LX
plsrca.res$ly / rplsres_all$lambda_0$LY

plsrca.res$fi.matrix / rplsres_all$lambda_0$fi
plsrca.res$fj.matrix / rplsres_all$lambda_0$fj
plsrca.res$p.matrix / rplsres_all$lambda_0$p
plsrca.res$q.matrix / rplsres_all$lambda_0$q
plsrca.res$pred.u.matrix / rplsres_all$lambda_0$pred.U
(plsrca.res$Zy.orig - plsrca.res$Y.rec) / rplsres_all$lambda_0$Y.resid
plsrca.res$Y.rec / rplsres_all$lambda_0$Y.rec

plsrca.res$X.recs / rplsres_all$lambda_0$X.recs
plsrca.res$X.resids / rplsres_all$lambda_0$X.resids

plsrca.res$Y.recs / rplsres_all$lambda_0$Y.recs
plsrca.res$Y.resids / rplsres_all$lambda_0$Y.resids


plsrca.res$r2.x.cumulative / rplsres_all$lambda_0$r2.x.cumulative
plsrca.res$r2.y.cumulative / rplsres_all$lambda_0$r2.y.cumulative
plsrca.res$r2.y / rplsres_all$lambda_0$r2.y
plsrca.res$r2.x / rplsres_all$lambda_0$r2.x



stupid.colors <- colorRampPalette(c("black","red"))(length(lambdas))
i <- 1
for(lambda in lambdas){
  
  takane.X <- X.res$Zx * sum(X)
  # takane.X <- diag(rowSums(X)) %*% (diag(1/X.res$m) %*% X.res$Zx)
  ### ORIGINALS... must reference these first.
  # takane.XM <- diag(rowSums(X)) + (lambda * matrix.generalized.inverse( tcrossprod(takane.X) ))
  # takane.XW <- diag(colSums(X)) + (lambda * ( t(takane.X) %*% matrix.generalized.inverse(tcrossprod(takane.X)) %*% takane.X ) )
  takane.XM <- diag(rowSums(X)) + (lambda * diag(nrow(X)))
  takane.XW <- diag(colSums(X)) + (lambda * diag(ncol(X)))
  # takane.XM <- rowSums(X) + (lambda * 1)
  # takane.XW <- colSums(X) + (lambda * 1)
  
  # print(diag(takane.XM) / takane.XM2)
  # print(diag(takane.XW) / takane.XW2)
  
  takane.Y <- Y.res$Zx * sum(Y)
  # takane.Y <- diag(rowSums(Y)) %*% (diag(1/Y.res$m) %*% Y.res$Zx)
  ### ORIGINALS... must reference these first.
  # takane.YM <- diag(rowSums(Y)) + (lambda * matrix.generalized.inverse( tcrossprod(takane.Y) ))
  # takane.YW <- diag(colSums(Y)) + (lambda * ( t(takane.Y) %*% matrix.generalized.inverse(tcrossprod(takane.Y)) %*% takane.Y ) )
  takane.YM <- diag(rowSums(Y)) + (lambda * diag(nrow(Y)))
  takane.YW <- diag(colSums(Y)) + (lambda * diag(ncol(Y)))
  # takane.YM <- rowSums(Y) + (lambda * 1)
  # takane.YW <- colSums(Y) + (lambda * 1)
  
  
  # print(diag(takane.YM) / takane.YM2)
  # print(diag(takane.YW) / takane.YW2)
  # 
  
  # print((takane.XM %^% (-1/2) %*% takane.X) / sweep(takane.X,1,sqrt(takane.XM2),'/') )
  # print(diag((takane.XW %^% (-1))) / (1/takane.XW2))
  
  # print((takane.YM %^% (-1/2)) / sweep(takane.Y,1,sqrt(1/takane.YM2),'*') )
  # print(diag((takane.YW %^% (-1))) / (1/takane.YW2))
  
  #pause()
  rplsres <- rplsca(X,Y,pls.type="reg",lambdas=lambda,ridge.type="fast",comps=0)
  takane.res <- gplsreg( ((takane.XM %^% (-1/2) %*% takane.X)  ), ( (takane.YM %^% (-1/2) %*% takane.Y)  ), takane.XW %^% (-1), takane.YW %^% (-1), comps = 0)
  
  plot.this <- takane.res$fj
  if(lambda==0){
    plot(plot.this,pch=20,col="black",xlim=c(-max(abs(plot.this[,1])*1.25),max(abs(plot.this[,1])*1.25)),ylim=c(-max(abs(plot.this[,2])*1.25),max(abs(plot.this[,2])*1.25)))
    #plot( (takane.res$d^2) / sum(takane.res$d^2),type="l",col=stupid.colors[i], ylim=c(0,1))
  }else{
    points(plot.this,pch=18,col=stupid.colors[i])
    #points( (takane.res$d^2) / sum(takane.res$d^2),type="l",col=stupid.colors[i])
  }
  # print(sum(takane.res$d^2) / sum(plsrca.res$d^2))
  
  # print(sum(rplsres[[1]]$d^2) / sum(takane.res$d^2))
  # print( head(rplsres[[1]]$fi / takane.res$fi) )
  # print( head(rplsres[[1]]$fj / takane.res$fj) )
  # pause()
  i <- i + 1
}



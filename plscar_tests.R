library(GSVD)
source("gplsreg_core.R")

## need to read in various types of data, call in from the PLSCAR paper and replce X & Y here.
library(ExPosition)


# data("TRAITS")
# data("SNPS")
# X <- make.data.nominal(TRAITS)
# Y <- make.data.nominal(SNPS)

 

load("/UTDallas/S17/Papers/PLSRCA/ADNI_DATA/TADPOLE.fin2.rda")
load("/UTDallas/S17/Papers/PLSRCA/ADNI_DATA/genetic.data.rda")
load("/UTDallas/S17/Papers/PLSRCA/ADNI_DATA/TADPOLE_Dict.fin2.rda")
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


CTX.uptake <- TADPOLE.fin[,grep("CTX",colnames(TADPOLE.fin))]
colnames(CTX.uptake) <- gsub("_UCBERKELEYAV45_10_17_16","",gsub("CTX_","",colnames(CTX.uptake)))

drop.rows <- which(rowSums(is.na(CTX.uptake))==70)
drop.cols <- grep("UNKNOWN",colnames(CTX.uptake))
CTX.uptake_drop <- CTX.uptake[-c(drop.rows),-c(drop.cols)]
CTX.uptake_drop <- apply(CTX.uptake_drop,2,function(x){ x[which(is.na(x))] <- mean(x,na.rm=T); x })

genetic.nom_drop <- genetic.nom[rownames(CTX.uptake_drop),]


X <- CTX.uptake_drop
Y <- genetic.nom_drop



# in this particular case it is assumed that X & Y will go through their ca preprocs on their own.
  Ox <- X/sum(X)
  mx <- rowSums(Ox)
  wx <- colSums(Ox)
  Ex <- mx %o% wx
  Zx <- Ox - Ex
  weightedZx <- sweep(sweep(Zx, 1, sqrt(mx), "/"), 2, sqrt(wx), "/")
  
  Oy <- Y/sum(Y)
  my <- rowSums(Oy)
  wy <- colSums(Oy)
  Ey <- my %o% wy
  Zy <- Oy - Ey
  weightedZy <- sweep(sweep(Zy, 1, sqrt(my), "/"), 2, sqrt(wy), "/")
  
  
  
  ## I can do this one in two ways: rely on the SVD optimization (V1) and the GSVD (V2) optimization
  ## I will actually need both for tests...
  
  start_1 <- Sys.time()
  gplsreg_results_V1 <- gpls_reg(X = weightedZx, Y = weightedZy)
  end_1 <- Sys.time()
  
  crossprod(gplsreg_results_V1$LX)
  diag(t(gplsreg_results_V1$LX) %*% gplsreg_results_V1$LY) / gplsreg_results_V1$Deltas
  crossprod(gplsreg_results_V1$Tmat)
  
  start_2 <- Sys.time()
  gplsreg_results_V2 <- gpls_reg(X = Zx, Y = Zy,
                                XLW = diag(1/mx), XRW = diag(1/wx),
                                YLW = diag(1/my), YRW = diag(1/wy)
                                )
  end_2 <- Sys.time()
  
  crossprod(gplsreg_results_V2$LX)
  diag(t(gplsreg_results_V2$LX) %*% gplsreg_results_V2$LY) / gplsreg_results_V1$Deltas
  crossprod(gplsreg_results_V2$Tmat)
  
  
    ### what is the relationship between these?
      ### it is something column-wise...
  gplsreg_results_V1$X_reconstructed[,,1] / gplsreg_results_V2$X_reconstructed[,,1]
  
    ## these are the same, this is where it comes from.
  (gplsreg_results_V1$pred_U / gplsreg_results_V2$pred_U)[,1]
  (gplsreg_results_V1$X_reconstructed[,,1] / gplsreg_results_V2$X_reconstructed[,,1])[1,]
    ## ok so I need to reconstruct pred_U from Tmat & the data to account for this?
      ## pred_U can be different I believe. Need to check with the nominals.
  
  
  ((gplsreg_results_V1$pred_U / gplsreg_results_V2$pred_U)[,1]) / ((gplsreg_results_V1$X_reconstructed[,,1] / gplsreg_results_V2$X_reconstructed[,,1])[1,])
  ((gplsreg_results_V1$pred_U / gplsreg_results_V2$pred_U)[,2]) / ((gplsreg_results_V1$X_reconstructed[,,2] / gplsreg_results_V2$X_reconstructed[,,2])[1,])
  ((gplsreg_results_V1$pred_U / gplsreg_results_V2$pred_U)[,3]) / ((gplsreg_results_V1$X_reconstructed[,,3] / gplsreg_results_V2$X_reconstructed[,,3])[1,])
  # pppaaaase()
  # 
  # gplsreg_results_V1$X_reconstructed[,,1] / ((diag(mx) %^% (-1)) %*% gplsreg_results_V2$X_reconstructed[,,1] %*% (diag(wx) %^% (-1)))
  
  gplsreg_results_V1$Y_reconstructed[,,1] / gplsreg_results_V2$Y_reconstructed[,,1]
  gplsreg_results_V1$Y_residuals[,,1] / gplsreg_results_V2$Y_residuals[,,1]
  
    ### the above tells me what's missing below.
  gplsreg_results_V1$X_residuals[,,1] / gplsreg_results_V2$X_residuals[,,1]
    ### close but not quite
    ((1/sqrt(mx)) %o% (1/sqrt(wx))) / (gplsreg_results_V1$X_residuals[,,1] / gplsreg_results_V2$X_residuals[,,1])
  
  
  # ((diag(1/mx) %^% (1/2)) %*% (Zx - gplsreg_results_V2$X_reconstructed[,,1]) %*% (diag(1/wx) %^% (1/2))) / gplsreg_results_V1$X_residuals[,,1]
  # (Zx - (diag(1/mx) %^% (-1/2)) %*% gplsreg_results_V2$X_reconstructed[,,1] %*% (diag(1/wx) %^% (-1/2)))  / gplsreg_results_V1$X_residuals[,,1]
  # 
  
  crossprod(gplsreg_results_V2$LX)
  diag(t(gplsreg_results_V2$LX) %*% gplsreg_results_V2$LY) / gplsreg_results_V2$Deltas
  crossprod(gplsreg_results_V2$Tmat)
  
  gplsreg_results_V1$LX / gplsreg_results_V2$LX
  gplsreg_results_V1$LY / gplsreg_results_V2$LY
  
  gplsreg_results_V1$Deltas / gplsreg_results_V2$Deltas
  gplsreg_results_V1$Betas / gplsreg_results_V2$Betas
  
  gplsreg_results_V1$U / gplsreg_results_V2$U
  gplsreg_results_V1$V / gplsreg_results_V2$V
  
  # gplsreg_results_V1$P / gplsreg_results_V2$P
  # gplsreg_results_V1$Q / gplsreg_results_V2$Q
  
  # gplsreg_results_V1$FI / gplsreg_results_V2$FI
  # gplsreg_results_V1$FJ / gplsreg_results_V2$FJ
  
  gplsreg_results_V1$Tmat / gplsreg_results_V2$Tmat
  # gplsreg_results_V1$pred_U / gplsreg_results_V2$pred_U
  
  gplsreg_results_V1$X_reconstructed / gplsreg_results_V2$X_reconstructed
  gplsreg_results_V1$X_residuals / gplsreg_results_V2$X_residuals
  
  gplsreg_results_V1$Y_reconstructed / gplsreg_results_V2$Y_reconstructed
  gplsreg_results_V1$Y_residuals / gplsreg_results_V2$Y_residuals
  
  gplsreg_results_V1$r2.x / gplsreg_results_V2$r2.x
  gplsreg_results_V1$r2.y / gplsreg_results_V2$r2.y
  
  
  
  ## now rebuild the Xs & Ys?
    ## start for now with just Ys.
  

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

#genetic.nom <- makeNominalData(genetic.data)
genetic.nom <- make_data_disjunctive(genetic.data)
genetic.nom.design <- make_data_disjunctive(as.matrix(attributes(genetic.nom)$variable.map))
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

DESIGN <- make_data_disjunctive(as.matrix(TADPOLE.fin$DX_bl))
  colnames(DESIGN) <- gsub("1\\.","",colnames(DESIGN))
W.DESIGN <- apply(DESIGN,2,function(x){x/sum(x)})
  rownames(W.DESIGN) <- rownames(DESIGN) <- rownames(TADPOLE.fin)
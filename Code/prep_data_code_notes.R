### some notes on the relative structure.

# ./plink --noweb --bfile "E:\Data\ADNI\Genetics\ADNI_GO_2_Forward_Bin" --bmerge "E:\Data\ADNI\Genetics\ADNI_GO2_GWAS_2nd_orig_BIN.bed" "E:\Data\ADNI\Genetics\ADNI_GO2_GWAS_2nd_orig_BIN.bim" "E:\Data\ADNI\Genetics\ADNI_GO2_GWAS_2nd_orig_BIN.fam" --make-bed --out "E:\Data\ADNI\Genetics\ADNI2MERGE"
# ./plink --noweb --bfile "E:\Data\ADNI\Genetics\ADNI2MERGE" --extract "E:\Data\Publications\PLSCAR\GETSNPS.txt" --maf 0.05 --mind 0.1 --geno 0.1 --recode --tab --out "E:\Data\Publications\PLSCAR\ADNIGo2_CANDIDATE_GENES"


## get reduced and preprocessed genetic data
library(data.table)
ped.file <- as.matrix(fread('E:/Data/Publications/PLSCAR/ADNIGo2_CANDIDATE_GENES.ped'))
ped.file <- gsub("[[:space:]]","",ped.file)
ped.ptids <- ped.file[,2]
ped.file <- ped.file[,-c(1:6)]
ped.file <- gsub("00",NA,ped.file)
ped.file <- cbind(ped.ptids,ped.file)
map.file <- as.matrix(fread('E:/Data/Publications/PLSCAR/ADNIGo2_CANDIDATE_GENES.map'))
colnames(ped.file) <- c("PTID",map.file[,2])

## clean up
genetic.data <- ped.file
rm(ped.file)
rm(map.file)


## read in TADPOLE data
TADPOLE_D1_D2_Dict <- read.csv('E:/Data/ADNI/Challenges/TADPOLE/TADPOLE_D1_D2_Dict.csv',header=T,stringsAsFactors = F)
TADPOLE_D1_D2 <- read.csv('E:/Data/ADNI/Challenges/TADPOLE/TADPOLE_D1_D2.csv',header=T,stringsAsFactors = F)

## match the IDs
overlapping.ids <- intersect(genetic.data[,1],TADPOLE_D1_D2$PTID)

## some clean up
ADNIGO.init.locs <- which( TADPOLE_D1_D2$ORIGPROT=="ADNIGO" & TADPOLE_D1_D2$VISCODE=="bl" & (TADPOLE_D1_D2$PTID %in% overlapping.ids)  )
ADNI2.init.locs <- which( TADPOLE_D1_D2$ORIGPROT=="ADNI2" & TADPOLE_D1_D2$VISCODE=="bl" & (TADPOLE_D1_D2$PTID %in% overlapping.ids)  )
TADPOLE_AGO2_BL <- TADPOLE_D1_D2[c(ADNIGO.init.locs,ADNI2.init.locs),]
TADPOLE_AGO2_BL <- TADPOLE_AGO2_BL[-c(which(colSums(is.na(TADPOLE_AGO2_BL))==nrow(TADPOLE_AGO2_BL)))]
TADPOLE_AGO2_BL <- TADPOLE_AGO2_BL[,-c(which((colSums(is.na(TADPOLE_AGO2_BL)) / nrow(TADPOLE_AGO2_BL)) > .1))]

## extract variables we want for analyses
AV45.vars <- TADPOLE_D1_D2_Dict[which(TADPOLE_D1_D2_Dict$TBLNAME=="UCBERKELEYAV45"),"FLDNAME"]
STRUCT.MRI.vars <- c("Ventricles_bl","Hippocampus_bl","WholeBrain_bl","ICV_bl")
KEY.DEMOG.vars <- c("RID","PTID","VISCODE","SITE","COLPROT","ORIGPROT","EXAMDATE","DX_bl","DX","DXCHANGE","AGE","PTGENDER","PTEDUCAT","PTETHCAT","PTRACCAT","APOE4")
CLIN.vars <- c("CDRSB","ADAS11","ADAS13","MMSE","MOCA",paste0(c("CDRSB","ADAS11","ADAS13","MMSE","MOCA"),"_bl"))
BROAD.NIMG.vars <- c("FDG","AV45",paste0(c("FDG","AV45"),"_bl"))
CTX.uptake.vars <- AV45.vars[which(grepl("CTX",AV45.vars) & (!grepl("SIZE",AV45.vars)))]

fin.vars <- c(KEY.DEMOG.vars,CLIN.vars,BROAD.NIMG.vars,CTX.uptake.vars,STRUCT.MRI.vars)

## prep final sets for saving
TADPOLE.fin <- TADPOLE_AGO2_BL[,fin.vars]
rownames(TADPOLE.fin) <- TADPOLE.fin$PTID
rownames(genetic.data) <- genetic.data[,"PTID"]
genetic.data <- genetic.data[rownames(TADPOLE.fin),]
TADPOLE_Dict.fin <- TADPOLE_D1_D2_Dict[which(TADPOLE_D1_D2_Dict$FLDNAME %in% colnames(TADPOLE.fin)),]

## save out
save(TADPOLE.fin,file="E:/Data/Publications/PLSCAR/TADPOLE.fin_2020DEC17.rda")
save(genetic.data,file="E:/Data/Publications/PLSCAR/genetic.data_Data_2020DEC17.rda")
save(TADPOLE_Dict.fin,file="E:/Data/Publications/PLSCAR/TADPOLE_Dict.fin_Data_2020DEC17.rda")


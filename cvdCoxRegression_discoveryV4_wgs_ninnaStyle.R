require(parallel)
require(data.table)
library(survival)
library(ISwR)
library(RColorBrewer)
library(RSQLite)
sqlite    <- dbDriver("SQLite")
shinyDB <- dbConnect(sqlite,"shinyKH.db")

nSites <- 10265398
Nind<-4629


## to get proper colours, does not import from .Rprofile for some reason
colorblind6b<-function(){
  grDevices::palette(c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
  grDevices::dev.off()
}
colorblind6b()


## getGeno function from /pontusData/greenland/shiny/greenlandWGSimpu2018e/functions.R
getGeno<-function(minSNP,numSNP,Nind,theFile="/home/albrecht/public/albrecht/shiny/greenland/data/geno.bin"){

  con<-file(theFile,"rb")
  seek(con,where=Nind*2*(minSNP-1))
  geno2<-readBin(con,"integer",Nind*numSNP,size=2)
  close(con)
  dim(geno2)<-c(Nind,length(geno2)/Nind)
  geno2[geno2==-1]<-NA

  geno2

}


dataFolder <- "/pontusData/greenland/shinyData/datagreenlandWGSimpu2018b/"

## indis ordered same way as
## /emc/emil/riccoStuff/gemmaTotal/data/wgsImpuDataV2/finalOnlyImpu.impute2_info080_wgsMergedV2_maf0005.fam
genoFile <- paste0( dataFolder, "geno.bin" )

nn<-load("/pontusData/greenland/shinyData/datagreenlandWGSimpu2018b/infoV3_RNA.Rdata")

fam<-fread2("/emc/emil/riccoStuff/gemmaTotal/data/wgsImpuDataV2/finalOnlyImpu.impute2_info080_wgsMergedV2_maf0005.fam",h=F)

## check that they are the same order
table(fam$V1==phe$particid)

## CVD registry data merged with data with PCs,sex,age and LDLR variant
k<-load("/home/emil/OMNI5/phenotypes/registryCVDdata270619includeDeathUntil2017withBirthyearAllAndBBHninnaStyle.Rdata")

bim<-fread2("/emc/emil/riccoStuff/gemmaTotal/data/wgsImpuDataV2/finalOnlyImpu.impute2_info080_wgsMergedV2_maf0005.bim",h=F)

int<-intersect(phe$particid,merged$particid)

merged<-merged[ order(match(merged$particid,phe$particid)),]


table(merged[ merged$particid%in%int,"particid"]==phe[ phe$particid%in%int,"particid"])


keep<-phe$particid%in%merged$particid

merged2<-merged[ merged$particid%in%int,]

table(merged2$particid==phe[ phe$particid%in%int,"particid"])

l1<-list(); l2<-list(); l3<-list(); l4<-list(); l5<-list(); l6<-list(); l7<-list()
r1<-list(); r2<-list(); r3<-list(); r4<-list(); r5<-list(); r6<-list(); r7<-list()


analyseCVD<-function(i,l1,l2,l3,l4,l5,l6,l7,r1,r2,r3,r4,r5,r6,r7,bim){

    if(i%%1000==0){
        print(i)
    }

    g<-getGeno( i, 1, Nind, genoFile )[keep]

    if(mean(as.numeric(g),na.rm = T)/2<0.05 | mean(as.numeric(g),na.rm = T)/2>0.95 | all(is.na(g)) | sum(is.na(g))/length(g)>0.01 ){

        l1[[bim[i,2]]]<-c(NA,NA,NA,NA)
        l2[[bim[i,2]]]<-c(NA,NA,NA,NA)
        l3[[bim[i,2]]]<-c(NA,NA,NA,NA)
        l4[[bim[i,2]]]<-c(NA,NA,NA,NA)
        l5[[bim[i,2]]]<-c(NA,NA,NA,NA)
        l6[[bim[i,2]]]<-c(NA,NA,NA,NA)
        l7[[bim[i,2]]]<-c(NA,NA,NA,NA)

        r1[[bim[i,2]]]<-c(NA,NA,NA,NA)
        r2[[bim[i,2]]]<-c(NA,NA,NA,NA)
        r3[[bim[i,2]]]<-c(NA,NA,NA,NA)
        r4[[bim[i,2]]]<-c(NA,NA,NA,NA)
        r5[[bim[i,2]]]<-c(NA,NA,NA,NA)
        r6[[bim[i,2]]]<-c(NA,NA,NA,NA)
        r7[[bim[i,2]]]<-c(NA,NA,NA,NA)

    } else{

        coxAll2<-coxph(formula = Surv(yearsToEvent,anyevent2017==1)~as.numeric(g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        l1[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventIHD,ihdevent2017==1)~as.numeric(g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        l2[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventHF,hfevent2017==1)~as.numeric(g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        l3[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventVASCUL,vasculevent2017==1)~as.numeric(g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        l4[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventKORONAR,koronarevent==1)~as.numeric(g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        l5[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventSTROKE,strokeevent2017==1)~as.numeric(g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        l7[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])

        ################### recessive

        coxAll2<-coxph(formula = Surv(yearsToEvent,anyevent2017==1)~as.numeric(2==g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        r1[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventIHD,ihdevent2017==1)~as.numeric(2==g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        r2[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventHF,hfevent2017==1)~as.numeric(2==g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        r3[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventVASCUL,vasculevent2017==1)~as.numeric(2==g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        r4[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventKORONAR,koronarevent2017==1)~as.numeric(2==g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        r5[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


        coxAll2<-coxph(formula = Surv(yearsToEventSTROKE,strokeevent2017==1)~as.numeric(2==g)+ALDER+KON+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=merged2)
        s<-summary(coxAll2)
        r7[[paste0(bim[i,2])]]<-c(s$coefficients[1,5],s$coefficients[1,2],bim[i,1],bim[i,4])


    }

    return(list(l1=l1,l2=l2,l3=l3,l4=l4,l5=l5,l7=l7,r1=r1,r2=r2,r3=r3,r4=r4,r5=r5,r7=r7))

}

m<-mclapply(1:nSites,function(x) analyseCVD(i=x,l1,l2,l3,l4,l5,l6,l7,r1,r2,r3,r4,r5,r6,r7,bim),mc.cores = 20)


save(m,file = paste0("cvd/coxRegCVD_discoveryAddRecTestV4_wgs_ninnaStyle_chrAll.Rdata") )
##save(l1,l2,l3,l4,l5,l7,file = paste0("coxRegCVD_discoveryAddV3_wgs_chrAll.Rdata") )
##save(r1,r2,r3,r4,r5,r7,file = paste0("coxRegCVD_discoveryRecV3_wgs_chrAll.Rdata") )




bitmap("cvd/coxRegCVD_discoveryV4addV2_wgs_ninnaStyle.png",res=200,width = 10)
par(mfrow=c(2,4));

n1<-sum(!is.na(lapply(l1,function(x) x[1])))
qqpemil(as.numeric(lapply(l1,function(x) x[1])),main="any CVD event");

n2<-sum(!is.na(lapply(l2,function(x) x[1])))
qqpemil(as.numeric(lapply(l2,function(x) x[1])),main="Ischemic heart disease");

n3<-sum(!is.na(lapply(l3,function(x) x[1])))
qqpemil(as.numeric(lapply(l3,function(x) x[1])),main="Heart failure");

n4<-sum(!is.na(lapply(l4,function(x) x[1])))
qqpemil(as.numeric(lapply(l4,function(x) x[1])),main="Peripheral artery disease");

n5<-sum(!is.na(lapply(l5,function(x) x[1])))
qqpemil(as.numeric(lapply(l5,function(x) x[1])),main="Coronary operations");

n6<-sum(!is.na(lapply(l6,function(x) x[1])))
qqpemil(as.numeric(lapply(l6,function(x) x[1])),main="Amputations");

n7<-sum(!is.na(lapply(l7,function(x) x[1])))
qqpemil(as.numeric(lapply(l7,function(x) x[1])),main="Stroke");

dev.off()

bitmap("wgs/coxRegCVD_discoveryV4recV2_wgs_ninnaStyle.png",res=200,width = 10)
par(mfrow=c(2,4));

n1<-sum(!is.na(lapply(l1,function(x) x[1])))
qqpemil(as.numeric(lapply(r1,function(x) x[1])),main="any CVD event");

n2<-sum(!is.na(lapply(l2,function(x) x[1])))
qqpemil(as.numeric(lapply(r2,function(x) x[1])),main="Ischemic heart disease");

n3<-sum(!is.na(lapply(l3,function(x) x[1])))
qqpemil(as.numeric(lapply(r3,function(x) x[1])),main="Heart failure");

n4<-sum(!is.na(lapply(l4,function(x) x[1])))
qqpemil(as.numeric(lapply(r4,function(x) x[1])),main="Peripheral artery disease");

n5<-sum(!is.na(lapply(l5,function(x) x[1])))
qqpemil(as.numeric(lapply(r5,function(x) x[1])),main="Coronary operations");

n6<-sum(!is.na(lapply(l6,function(x) x[1])))
qqpemil(as.numeric(lapply(r6,function(x) x[1])),main="Amputations");

n7<-sum(!is.na(lapply(l7,function(x) x[1])))
qqpemil(as.numeric(lapply(r7,function(x) x[1])),main="Stroke");

dev.off()

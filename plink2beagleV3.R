args<-commandArgs(trailingOnly=T)

if(length(args)==0){
  print("Converts a binary plink file (.bed,.fam,.bim) to a (TAB-delimeted) beagle file with genotype likelihoods - for a specified depth")
  print("1. plink File, 2. depth, 3. how many individuals in each depth group seperated by ',' 4. .frq file from plink")
  print("It is assumed that the first allele of the .bim file is coded as 0 and the second allele coded as 1 in the .bed file")
  print("The marker ID of the beagle file will be chr_pos from the .bim file")
  q()
}

## SEE PLINK DOCUMENTATION ABOUT CODING OF .bed FILE
## https://www.cog-genomics.org/plink/1.9/formats#bed

plinkFile<-args[1]
depth<-args[2]

pl<-plinkV2(plinkFile)
##pl<-plinkV2("/home/emil/OMNI5/ADMIXED/OMNI/datPlus_removeRelatedHet_maf005_geno005_thin01")

if(length(args)>2){
    Nind<-args[3]
    vecNind<-as.numeric(unlist(strsplit(Nind,",")))
    if(nrow(pl$geno)!=sum(vecNind)){
        print("Number of individuals specified not the same as number of individuals in plink file!")
        stop() 
    }
}

if(length(args)>3){    
    frqFile<-args[4]
    frq<-read.table(frqFile,as.is=T,h=T)
}

if(length(args)>2){

    depth<-as.numeric(unlist(strsplit(depth,",")))
    depth<-rep(c(depth),times=c(vecNind))
}

depth<-as.numeric(depth)

##############

#x=genotype
#d=dybde (avg)
#e=error rate
#norm=Normalize to sum to one
getLikes<-function(x,d=2,e=0.01,norm=FALSE){
  n<-length(x)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[x+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm){
    res<-t(t(res)/colSums(res))
  }
  res[is.na(res)]<-1/3
  return(res)
}





## swap around genotypes so 2 is minor/minor (counting minors)
k<-t(apply(pl$geno,2,function(y) getLikes(x=(2-y),d=depth,norm=TRUE)))

## with MAF prior


if(length(args)>3){
    
    kMaf<-t(sapply(1:nrow(k),function(x) k[x,]*rep(c((1-frq[x,"MAF"])**2,2*(1-frq[x,"MAF"])*frq[x,"MAF"],frq[x,"MAF"]**2),ncol(k)/3)))
            
    norm<-t(apply(kMaf,1,function(x) sapply(seq(from=3,to=3*ncol(kMaf),by=3),function(y) x[y-2]+x[y-1]+x[y])))
    
    kMafNorm<-kMaf/norm

    kMaf2<-cbind(paste0(pl$bim$V1,"_",pl$bim$V4),pl$bim$V6,pl$bim$V5,kMafNorm)

    colnames(kMaf2)<-c("marker","allele1","allele2",paste0(rep("Ind",3),sapply(1:nrow(pl$geno),function(x) rep(x,3))))
    
    write.table(kMaf2,paste0(plinkFile,"_depth",depth,".mafprior.beagle"),col=T,row=F,qu=F,sep = "\t")

} else if(length(args)>2){

    k2<-cbind(paste0(pl$bim$V1,"_",pl$bim$V4),pl$bim$V6,pl$bim$V5,k)

    colnames(k2)<-c("marker","allele1","allele2",paste0(rep("Ind",3),sapply(1:nrow(pl$geno),function(x) rep(x,3))))

    depth<-paste0(unique(depth),collapse="_")
    
    write.table(k2,paste0(plinkFile,"_depth",depth,".beagle"),col=T,row=F,qu=F,sep = "\t")
    
    ## with uniform prior
} else{

    k2<-cbind(paste0(pl$bim$V1,"_",pl$bim$V4),pl$bim$V6,pl$bim$V5,k)

    colnames(k2)<-c("marker","allele1","allele2",paste0(rep("Ind",3),sapply(1:nrow(pl$geno),function(x) rep(x,3))))
    
    write.table(k2,paste0(plinkFile,"_depth",depth,".beagle"),col=T,row=F,qu=F,sep = "\t")
}

##apply(k2[,4:ncol(k)],1,function(x)   )

## check that we get the same as genotypes
##site<-2
##g<-sapply(1:nrow(pl$geno),function(x) {print(x); which.max(c(as.numeric(k2[site,(3*x+1)]),as.numeric(k2[site,(3*x+2)]),
##                                                             as.numeric(k2[site,(3*x+3)])))})

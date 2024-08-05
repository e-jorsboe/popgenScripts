library(parallel)
require(data.table)

arguments<-commandArgs(trailingOnly=T)

bimFile1<-arguments[1]
bimFile2<-arguments[2]
cores<-as.numeric(arguments[3])


if(length(arguments)==0){
  print("1st argument is first bim file, 2nd argument is second bim file, finds overlapping sites chr_pos_A0_A1")
  q()
  
}

if(!length(arguments)>1 | !length(arguments)>2){
    print("has to give 2nd argument, second bim file")
    print("has to give 3rd argument, number of cores")
    q()
  
}


bim1<-fread(bimFile1,data.table=F,h=F)
bim2<-fread(bimFile2,data.table=F,h=F)

fun2<-function(x){return(paste(sort(c(a0[x],a1[x])),collapse="_"))}

a0<-as.vector(bim1$V5)
a1<-as.vector(bim1$V6)
alleles1<-parallel::mclapply(1:nrow(bim1),fun2,mc.cores=cores)
bim1$id<-paste0(bim1$V1,"_",bim1$V4,"_",alleles1)

a0<-as.vector(bim2$V5)
a1<-as.vector(bim2$V6)
alleles2<-parallel::mclapply(1:nrow(bim2),fun2,mc.cores=cores)
bim2$id<-paste0(bim2$V1,"_",bim2$V4,"_",alleles2)


print("based on chr_pos id")
print(length(intersect(paste0(bim1$V1,bim1$V4),paste0(bim2$V1,bim2$V4))))

print("based on chr_pos_A0_A1 id")
print(length(intersect(bim1$id,bim2$id)))

bimFileName1<-rev(unlist(strsplit(bimFile1,split = "/")))[1]
bimFileName2<-rev(unlist(strsplit(bimFile2,split = "/")))[1]


print("Removes duplicated IDs from both datasets, as these are essentially the same site")

bim1V2<-bim1[ !duplicated(bim1$id),]
bim2V2<-bim2[ !duplicated(bim2$id),]


int<-intersect(bim1V2$id,bim2$id)

bim1V3<-bim1V2[ bim1V2$id%in%int,]
bim2V3<-bim2V2[ bim2V2$id%in%int,]

bim2V3<-bim2V3[ order(match(bim2V3$id,bim1V3$id)),]

if(!all(bim1V3$id==bim2V3$id)){
    q("bim files not sorted the same way - something is wrong!")
}

print("After removing duplicate sites")
print(length(intersect(bim1V3$id,bim2V3$id)))

overlap<-cbind(bim1V3,bim2V3[,2])
colnames(overlap)<-c("CHR","ID1","MORGAN","POS","A0","A1","ID2")

write.table(overlap,paste0("overlap_",bimFileName1,bimFileName2,".bim"),col=F,row=F,quote = F)

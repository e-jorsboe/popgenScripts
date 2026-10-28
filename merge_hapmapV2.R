bimFile<-commandArgs(trailingOnly=T)[1]
## which chr
i<-as.numeric(commandArgs(trailingOnly=T)[2])
prefix<-commandArgs(trailingOnly=T)[3]

bim<-read.table(bimFile, head=F,as.is=T)

## finds overlaps and extract info on overlapped markers

## so those sites that are in both our data and the hapmap data are extraced



hap<-read.table(paste("/home/emil/OMNI5/genetic_map_b37/genetic_map_chr",i,"_combined_b37.txt",sep=""),head=T,as.is=T)

int<-intersect(hap$position,bim$V4)


if(all(bim[ bim$V4%in%int,"V4"]==hap[ hap$position%in%int,"position"])){

    bim2<-bim[ bim$V4%in%int,1:3]
    bim3<-cbind(bim[,4],0,0)
    bim3<-data.frame(bim3,stringsAsFactors=F)
    bim2[,1:3]<-hap[ hap$position%in%int,1:3]

    if(all(bim3[ bim3[,1]%in%int,1]==bim2[,1])){
        bim3[ bim3[,1]%in%int,2:3]<-bim2[,2:3]
    }
    colnames(bim3)<-colnames(hap)
}

write.table(bim3,file=paste(prefix,"_map_chr",i,"_b37.txt",sep=""),row.names=F,col.names = T,quote=F,sep=" ")
print(i)

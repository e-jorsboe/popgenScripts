arguments<-commandArgs(trailingOnly=T)

Qfile<-arguments[1]


if(length(arguments)==0){
  print("1st argument is .Q file, 2nd argument (OPTIONAL) is list with group names for each individual - MUST BE SORTED same way!")
  print("If second argument is supplied will generate admixture barplot grouped by provided grouping - sorted on Q column with highest mean in each group")
  q()
  
}

Q<-read.table(Qfile,as.is=T)

pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")

palette(pal(10))


if(length(arguments)>1){
  groupsFile<-arguments[2]
  groups<-read.table(groupsFile,as.is=T)
  ## want to see these values as names/labels rather than for instance integers
  groups$V1<-as.character(groups$V1)
  
  if(nrow(Q)!=nrow(groups)){
    print("2 Files must have same number of rows!!")
    q()
  }
  
  
  
}  

##find which one Q column has highest mean value
maxCol<-which.max(colMeans(Q))
if(length(arguments)>1){
  sorting<-order(groups$V1,Q[,maxCol])
  Q<-Q[ sorting,]
  groups<-data.frame(V1=groups[ sorting,],stringsAsFactors = F)
} else{
  
  Q<-Q[ order(Q[,maxCol]),]
}


tQ<-t(Q)

tmp<-tQ[1,]
tQ[1,]<-tQ[maxCol,]
tQ[maxCol,]<-tmp


if(length(arguments)>1){
  
  bitmap(paste(tail(unlist(strsplit(Qfile,"/")),n=1),"_grouped.png",sep=""),res=300,width=12)
  par(mar=c(3, 2.5, 2, 1))
  
  h<-barplot(tQ,col=1:nrow(tQ),border=NA,space=0,ylab="admixture proportion",main=paste("K=",nrow(tQ),sep=""),names.arg = rep("",ncol(tQ)))
  abline(v=h[cumsum(table(groups$V1))][-length(unique(groups$V1))]+0.5)
  text(tapply(h,groups$V1,mean),rep(-0.05,nrow(tQ)),unique(groups$V1),xpd=TRUE)
  
  dev.off()
  
} else{
  
  bitmap(paste(tail(unlist(strsplit(Qfile,"/")),n=1),".png",sep=""),res=300,width=12)
  par(mar=c(3, 2.5, 2, 1))
  
  h<-barplot(tQ,col=1:nrow(tQ),border=NA,space=0,ylab="admixture proportion",main=paste("K=",nrow(tQ),sep=""))
  
  dev.off()
  
}

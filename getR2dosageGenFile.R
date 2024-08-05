require(data.table)


args<-commandArgs(trailingOnly = T)

genfile<-args[1]
SNP<-args[2]

gen<-fread2(genfile,h=F)

rows<-(ncol(gen)-5)/3

a<-t(apply(gen[,6:ncol(gen)],1,function(x) rowSums(matrix(x*rep(0:2,times=rows),nrow=rows,ncol=3,byrow=T))))

r2<-apply(a,1, function(x) cor(x,a[which(gen$V2%in%SNP),],use="pairwise.complete.obs")**2)

df<-data.frame(rsID=gen$V2,r2=as.numeric(r2),stringsAsFactors=F)

df2<-df[ df$r2>0.25,]

write.table(df,paste0(genfile,".doseLD"),row.names=T,qu=F,col.names=F)

write.table(df2,paste0(genfile,"R2025.doseLD"),row.names=T,qu=F,col.names=F)

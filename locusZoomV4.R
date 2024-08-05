library(RColorBrewer)

l<-commandArgs(trailingOnly = TRUE)

pl  <- l[1]
inputChr <- as.numeric(l[2])
begin_pos <- as.numeric(l[3])
end_pos <- as.numeric(l[4])
asso_name <- l[5]
pheno <- l[6]
snp <- l[7]

##pl<-"/emc/emil/riccoStuff/gemmaTotal/data/locusZoom"
##inputChr<-11
##begin_pos <- 116000000 
##end_pos <- 118000000

##pheno<-"apoA1"
## "output/" added in bash script
##asso_name<-"output/run2017_v4_addtrans_ALL_C18_1_n9c.assoc.txt"
##asso_name<-"output/run2017_v5_addtrans_ALL_apoA1.assoc.txt"
##snp<-"11:117063003-G-A"
##l<-list(); for(i in 1:7){l[[i]]<-i}

print(length(l))

if(length(l) > 6 & snp!="0"){
  run<-paste0(unlist(strsplit(basename(asso_name),"_"))[1:2],collapse = "_")
  print(run)
  model<-unlist(strsplit(basename(asso_name),"_"))[3]
  print(model)
  cohort<-unlist(strsplit(basename(asso_name),"_"))[4]
  print(cohort)
  pl2<-plinkV2("/emc/emil/riccoStuff/gemmaTotal/data/snp")
  cov<-read.table(paste0("/emc/emil/riccoStuff/gemmaTotal/",run,"/",pheno,"/",run,"_",model,"_",cohort,"_",pheno,".cov"),as.is=T)
  cov2<-cbind(cov,pl2$geno[,1])
  write.table(cov2,"/emc/emil/riccoStuff/gemmaTotal/data/snp.cov",col=F,row=F,qu=F)
  if(grepl("rec",model)){
    print(paste0("/emc/emil/riccoStuff/gemmaTotal/prog/gemma_recessive -bfile /emc/emil/riccoStuff/gemmaTotal/data/locusZoomConditional -k /emc/emil/riccoStuff/gemmaTotal/",sub("rectrans_","",sub("assoc","sXX",asso_name))," -lmm 4 -o conditional_",snp,"_",run,"_",model,"_",cohort,"_",pheno," -c /emc/emil/riccoStuff/gemmaTotal/data/snp.cov -maf 0 -miss 0.8"))
    system(paste0("/emc/emil/riccoStuff/gemmaTotal/prog/gemma_recessive -bfile /emc/emil/riccoStuff/gemmaTotal/data/locusZoomConditional -k /emc/emil/riccoStuff/gemmaTotal/",sub("rectrans_","",sub("assoc","sXX",asso_name))," -lmm 4 -o conditional_",snp,"_",run,"_",model,"_",cohort,"_",pheno," -c /emc/emil/riccoStuff/gemmaTotal/data/snp.cov -maf 0 -miss 0.8"))
  } else{
    print(paste0("/emc/emil/riccoStuff/gemmaTotal/prog/gemma.linux -bfile /emc/emil/riccoStuff/gemmaTotal/data/locusZoomConditional -k /emc/emil/riccoStuff/gemmaTotal/",sub("addtrans_","",sub("assoc","sXX",asso_name))," -lmm 4 -o conditional_",snp,"_",run,"_",model,"_",cohort,"_",pheno," -c /emc/emil/riccoStuff/gemmaTotal/data/snp.cov -maf 0 -miss 0.8"))
    system(paste0("/emc/emil/riccoStuff/gemmaTotal/prog/gemma.linux -bfile /emc/emil/riccoStuff/gemmaTotal/data/locusZoomConditional -k /emc/emil/riccoStuff/gemmaTotal/",sub("addtrans_","",sub("assoc","sXX",asso_name))," -lmm 4 -o conditional_",snp,"_",run,"_",model,"_",cohort,"_",pheno," -c /emc/emil/riccoStuff/gemmaTotal/data/snp.cov -maf 0 -miss 0.8"))
  }
}

asso<-read.table(asso_name,head=T,as.is=T,comment.char="")

################################################################################################################
## IF .ASSOC FILE DOESNOT HAVE "p_lrt" COLUMN MIGHT BE AN ISSUE - AS USED LATER IN SCRIPT FOR GETTING P VALUES
###################################################################################################################

## if .assoc only has p_score and not p_lrt column (assumed to be column with P values in this script)
if(any(grepl("p_score",colnames(asso))) & (all(!grepl("p_lrt",colnames(asso))) & all(!grepl("p_wald",colnames(asso))))){
  print("Only has p_score, will base P values on that!")
  colnames(asso)[which(colnames(asso)=="p_score")]<-"p_lrt"
}

if(length(l) > 6  & snp!="0"){
  the_snp<-asso[asso$rs==snp,]
  print(paste0("output/conditional_",snp,"_",run,"_",model,"_",cohort,"_",pheno,".assoc.txt"))
  assoCond<-read.table(paste0("output/conditional_",snp,"_",run,"_",model,"_",cohort,"_",pheno,".assoc.txt"),as.is=T,h=T)
  asso2<-rbind(assoCond,the_snp)  
  asso2<-asso2[order(asso2$ps),]
} else{
  asso2<-asso  
}

plin<-plinkV2(pl)
geno<-plin$geno

bim<-plin$bim
n1<-bim[bim$V4>=begin_pos & bim$V4<=end_pos & bim$V1==inputChr,"V2"]

asso3<-asso2[asso2$rs %in% n1,]
n2<-n1 %in% asso3$rs
n3<-n1[n2]
geno2<-geno[,n3]

print("asso3:")
print(dim(asso3))

print("geno2:")
print(dim(geno2))

##    geno is a Nind x Nsnp matrix (\in 0,1,2,NA)
### Take like from plink geno

##    pval is a Nsnp vector of pvalues
### THere needs to be a pval for each SNP

##    pos is a Nsnp vector of positions
### same as pval

##    chr is a Nsnp vector of chromosomes
### same as pval

##    w indicate which SNPs is the SNP in focus (e.g. 5)
### the number of the SNP in focus

##    rs is the name of the SNP or the title of the plot
    
locusZoom<-function(geno,pval,pos,chr,w,rs="my title",con_snp=NULL){
    xforrs = 0.03
    regsz = 1.2
    width=22
    rsqplus = 0.045
    rightylabxplus=0.05
    xmargin = 0.005
    cexaxis=1.5
    cexlab=1.5
    pval[is.na(pval)]<-1  
    N <- length(pval)
    ld <- Relate:::ld.snp3(geno+1,back=N-1)
    rsqwithn = c(rev(sapply(1:(w-1),function(i){ld$rmisc[i,w-i]})),1,ld$rmisc[1:(N-w),w])
    palette(brewer.pal(10, "RdBu"))[1:10]
    blue="dodgerblue4"
    cols = c(10,8,6,4,2,2)
    getcol <- function(x){cuts = c(.2,.4,.6,.8,1);cols[6-sum(as.numeric(x)<cuts)]}
    min.pos <- min(pos)/1e6
    max.pos <- max(pos)/1e6
    rec = read.table(paste("/home/emil/OMNI5/genetic_map_b37/forLocusZoom/genetic_map_GRCh37_chr",chr[w],".txt",sep=""),header=T)
    relrec = rec[(rec[,2]/1e6)>min.pos&(rec[,2]/1e6)<max.pos,] 
    adj = 0
    nf <- layout(matrix(c(1,2), 2, 1, byrow=TRUE),heights=c(5.5,2.65),widths=c(8))
    par(mar=c(0.2, 6.2, 4, 6.3) + 0.1)
    plot(relrec[,"Position.bp."]/1e6,relrec[,"Rate.cM.Mb."],type="l",ylim=c(0,105),xaxs="i", xlim=c(min.pos-adj,max.pos+adj),axes=F,col=blue,xlab="",ylab="",lwd=3.5,cex.axis=cexaxis)
    axis(4, labels=seq(0,100,20),at=seq(0,100,20),col=blue,las=1,col.axis=blue,cex.axis=cexaxis)
    ##  text(par("usr")[2]+rightylabxplus,45,xpd=TRUE,srt=-90,labels=c(expression(paste("Recombination rate (cM ", Mb^-1,")"))),col=blue,cex=cexlab)
    mtext(c(expression(paste("Recombination rate (cM ", Mb^-1,")"))),4,4,col=blue,cex=cexlab)
    par(new=T)
    bgcols=sapply(rsqwithn,getcol)
    bgcols[w]=1
    ylim=c(0,max(10.5,max(-log10(pval))))
    plot(pos/1e6,-log10(pval),col="black",bg=bgcols,xaxs="i",xlim=c(min.pos-adj,max.pos+adj),xaxt="n",xlab="",ylab="",ylim=ylim,pch=21,cex=1,las=1,cex.axis=cexaxis,main=rs)
    ## add circle around SNP conditioned on
    if(!is.null(con_snp)){
      points(asso3[ asso3$rs==con_snp,"ps"]/1e6, -log10(asso3[ asso3$rs==con_snp,"p_lrt"]),pch=21,cex=2) 
    }
    axis(2,ylim=c(0,15),col="black",label=FALSE,cex.axis=cexaxis)
    mtext(2,text=expression(paste("-l",og[10],"[",italic(P),"]")),line=3,cex=cexlab)
    xstart = max.pos - (max.pos-min.pos)*0.15
    scalewidth=(max.pos-min.pos)*0.0125
    xx = c(xstart,xstart,xstart+scalewidth,xstart+scalewidth)
    ##  legend("topleft",legend=c(rs),pt.bg=1,col="black",cex=1.6,pch=21)
    ybot=ylim[2]*0.65
    ytop=ylim[2]*0.95
    ysz =ytop-ybot
    cuts = c(.2,.4,.6,.8,1)
    txtcuts = c( "0.2","0.4","0.6","0.8","1.0")
    txtposs = c(0.2,0.4,0.6,0.8,1.00)
    scalexplus=(max.pos-min.pos)*0.06
    polygon(xx+scalexplus,c(ybot,ybot+cuts[1]*ysz,ybot+cuts[1]*ysz,ybot),col=cols[1])
    for(i in 2:length(cols)){
        polygon(xx+scalexplus,c(ybot+cuts[(i-1)]*ysz,ybot+cuts[i]*ysz,ybot+cuts[i]*ysz,ybot+ysz*cuts[(i-1)]),col=cols[i])
    }
    scalenumplus =(max.pos-min.pos)*0.09
    ## text((xx[1]+0.4+xx[3])/2+rsqplus,ytop+0.35,expression(italic(r)^2),cex=cexaxis)
    text(xx[3]+scalenumplus,ytop+0.35,expression(italic(r)^2),cex=cexaxis)
    text(rep(xx[3]+scalenumplus,2),txtposs[3:4]*ysz+ybot-0.07,txtcuts[3:4],cex=cexlab-0.1)
    text(rep(xx[3]+scalenumplus,1),txtposs[2]*ysz+ybot-0.07,txtcuts[2],cex=cexlab-0.1)
    text(rep(xx[3]+scalenumplus,1),txtposs[1]*ysz+ybot-0.07,txtcuts[1],cex=cexlab-0.1)   
    dat<-read.table("/home/emil/OMNI5/refGeneHG19.gz",as.is=T,head=T,comment.char="")
    xx2 = dat[dat[,"chrom"]==paste("chr",chr[w],sep="") & dat[,"cdsStart"]<max.pos*1e6 & dat[,"cdsEnd"] >min.pos*1e6,]
    start = xx2$txStart
    end   = xx2$txEnd
    nams  = xx2$name2
    cnts  = xx2$exonCount
    abline(v=pos,lty=2)
    par(mar=c(5.2, 6.2, -0.1, 6.3) +0.1)
    plot(c(0,0),c(0,0),type="n",xlim=c(min.pos-adj,max.pos+adj),ylim=c(-0.8,0.1),xlab="",xaxs="i",yaxt='n',ylab="",main="",cex.lab=2.6,cex.axis=cexaxis,tck=-0.05)
    mtext(1,text=paste("Position on chromosome ",chr," (Mb)",sep=""),line=3,cex=cexlab)
    abline(v=pos,lty=2)
    ord <- order(start)
    start    <- start[ord]
    end      <- end[ord]
    exoncnts <- cnts[ord]
    nams     <- nams[ord]
    keep <- !duplicated(nams)
    start    <- start[keep]
    end      <- end[keep]
    exoncnts <- cnts[keep]
    nams     <- nams[keep]
    ord <- ord[keep]
    he       <- rep(c(0,-0.18,-0.36,-0.54,-0.72),100)[1:length(nams)]-0.05
if(length(start)>0){
    segments(start/1e6, he, end/1e6, he)
    keep = !duplicated(nams)
    sapply(1:sum(keep),function(x){text((end[keep][x]+start[keep][x])/2e6,he[keep][x]+0.08,bquote(italic(.(nams[keep][x]))),cex=cexlab-0.6)})
    estart = as.numeric(unlist(sapply(xx2$exonStarts[ord],function(y){strsplit(y,",")[[1]]})))/1e6
    eend = as.numeric(unlist(sapply(xx2$exonEnds[ord],function(y){strsplit(y,",")[[1]]})))/1e6
    rect(estart,rep(he,xx2$exonCount[ord])-0.01,eend, rep(he,xx2$exonCount[ord])+0.01,col="black")
}    
}

if(length(l) > 6 & snp!="0"){
    png(paste("conditional_",snp,"_",pheno,"_",begin_pos,"_",end_pos,"_chr",inputChr,".png",sep=""),width=1000)
    locusZoom(geno2,asso3$p_lrt,asso3$ps,rep(inputChr,length(n3)),sort.list(asso3$p_lrt)[1],paste("conditional_",the_snp$rs,"_",pheno,"_",begin_pos,"_",end_pos,"_chr",inputChr,sep=""),con_snp=snp)
    dev.off()
} else{
    png(paste(pheno,"_",begin_pos,"_",end_pos,"_chr",inputChr,".png",sep=""),width=1000)
    locusZoom(geno2,asso3$p_lrt,asso3$ps,rep(inputChr,length(n3)),sort.list(asso3$p_lrt)[1],paste(pheno,"_",begin_pos,"_",end_pos,"_chr",inputChr,sep=""))
    dev.off() 
}

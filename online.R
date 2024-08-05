winSum<-function(x,step){
    cs<-c(0,cumsum(x))
    cs[-c(1:step)]-cs[-(length(cs) - 1:step+1)]
}

winMean<-function(x,step){
    cs<-c(0,cumsum(x))
    (cs[-c(1:step)]-cs[-(length(cs) - 1:step+1)])/step
}


winMedian<-function(x,step){
    if(step>length(x)){
        print("step bigger than input")
        return()
    }
    return(sapply(step:length(x), function(y) median(tail(head(x,n=y),n=step))))
}


latticeplot<-function(mat,col=rev(heat.colors(50))[-1],at=0:10/10,...)
  lattice::levelplot(mat,scales=list(x=list(rot=45)),col.regions=col,colorkey=list(at=at,label=list(at=at,lab=at)),at=at,...)

info<-function(){
cat("xpd=T\n")
cat("par(mar=c(5.1,4.1,4.1,2.1))\n")
}
pdff<-function(...)
pdf("~/public/albrecht/temp.pdf",...)

norm<-function(x)
  x/sum(x)

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

brewer<-function(x="Dark2"){
    
pdf("/dev/null")
  library(RColorBrewer)
  palette(brewer.pal(8, x))
dev.off()

}

color.palette <- function(steps, n.steps.between=NULL, ...){

    if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
    if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")

    fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)

    for(i in which(n.steps.between>0)){
        col.start=RGB[,fill.steps[i]]
        col.end=RGB[,fill.steps[i+1]]
        for(j in seq(3)){
            vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]
            RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
        }
    }

    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}



myPals<-function(x=0){

    if(x==0)
        pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")
    else{
        
        require(wesanderson)
        
        aa<-c( "GrandBudapest","Moonrise1", "Royal1", "Moonrise2", "Cavalcanti", "Royal2","GrandBudapest2", "Moonrise3", "Chevalier" , "BottleRocket" ,"darjeeling", "darjeeling2")
        pal <- wes.palette(name = aa[x], type = "continuous")
    }
    
    pal
}

homo<-function()
  palette(c("mistyrose","lavender","lightyellow","lightblue","lightgreen","seashell","lightcyan"))

qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}

qqp<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,col=1,...){
  x<-x[!is.na(x)]
  if(!missing(maxLogP))
    x[x<10^-maxLogP]<-10^-maxLogP
  N<-length(x)
 chi1<-qchisq(1-x,1)
 x<-sort(x)
 lambda<-round(median(chi1)/qchisq(0.5,1),2)
  e<- -log((1:N-0.5)/N,10)
  if(add)
    points(e,-log(x,10),...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
    abline(0,1,col=2,lwd=2)
  }
  if(ci){
    c95<-qbeta(0.95,1:N,N-(1:N)+1)
    c05<-qbeta(0.05,1:N,N-(1:N)+1) 
    lines(e,-log(c95,10))
    lines(e,-log(c05,10))
  }   
}

qqpemil<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,...){
    x<-x[!is.na(x)]
    if(!missing(maxLogP)){
        x[x<10^-maxLogP]<-10^-maxLogP
    }
    N<-length(x)
    chi1<-qchisq(1-x,1)
    x<-sort(x)
    lambda<-round(median(chi1)/qchisq(0.5,1),2)
    e<- -log((1:N-0.5)/N,10)
    if(add){
        points(e,-log(x,10),...)
    } else{
        plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
        abline(0,1,col=2,lwd=2)
    }
    
    ##title(paste("lambda=",lambda), cex=1.5)
    mtext(paste0("lamdbda = ",lambda))
    
    if(ci){
        c95<-qbeta(0.95,1:N,N-(1:N)+1)
        c05<-qbeta(0.05,1:N,N-(1:N)+1) 
        lines(e,-log(c95,10))
        lines(e,-log(c05,10))
    }
}



bases<-c("A","C","G","T")
GENO<-c("AA","AC","AG","AT","CC","CG","CT","GG","GT","TT")
GENO2<-c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
GENOrev<-c("TT","GT","CT","AT","GG","CG","AG","CC","AC","AA")
GENO2rev<-rev(GENO2)

wiki<-function(x,file,...){
  
  if(class(x)=="data.frame")
    x<-as.matrix(x)
  if(class(x)=="matrix")
    wiki.matrix(x,file,...)
  


}

wiki.matrix<-function(x,file,class="wikitable",border="1",...){
write(paste("{| class=\"",class,"\" border=\"",border,"\"",sep=""),file=file)
temp<-NULL

if(!is.null(rownames(x)))
  temp<-" "
if(!is.null(colnames(x)))
  writeWikiRow(colnames(x),file=file,type="!",row=temp)

for(tal in 1:nrow(x))
  writeWikiRow(x[tal,],file=file,row=rownames(x)[tal])


write("|}",file=file,append=TRUE)
}

writeWikiRow<-function(x,file,type="|",row=NULL){
  write("|-",file=file,append=TRUE)
  if(!is.null(row))
  write(paste("!",row),file=file,append=TRUE)
  for(tal in 1:length(x))
  write(paste(type,x[tal]),file=file,append=TRUE)
}
wikii<-function(x,...)
wiki(x,file="~/public/albrecht/temp.txt",...)


#################################################
swichMe<-function(w,x){
  n<-length(w)
  if(n==nrow(x)){
    x<-t(x)
    warning("transposing x (slower)")
  }
  if(n!=ncol(x)){
    stop("wrong Dim")
  } 
  x[w+(0:(n-1))*nrow(x)]
}

qtrans<-function(x){
    k<-!is.na(x)
    ran<-rank(x[k])
    y<-qnorm((1:sum(k)-0.5)/sum(k))
    x[k]<-y[ran]
    x
}

##################################

plink<-function(plinkFile){
    pl<-snpMatrix::read.plink(plinkFile)
    ind<-rownames(pl)
    snp<-colnames(pl)
    geno<-as.integer(as.integer(pl)-1)
    dim(geno)<-dim(pl)
    geno[geno==-1]<-NA
    rownames(geno)<-ind
    colnames(geno)<-colnames(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    list(geno=geno,bim=bim,fam=fam,pl=pl)
}

## DO NOT USE - uses too much RAM
plinkV2<-function(plinkFile){
    pl<-snpStats::read.plink(plinkFile)
    pl2<-matrix(methods::as(pl$genotypes,"numeric"),nrow=nrow(pl$genotypes),ncol=ncol(pl$genotypes))
    colnames(pl2)<-colnames(pl$genotypes)
    
    ind<-rownames(pl2)
    snp<-colnames(pl2)
    ##    geno<-as.integer(as.integer(pl)-1)
    ##dim(geno)<-dim(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    rownames(pl2)<-fam$V2
    ## pl is a snpMatrix object might be different
    list(geno=pl2,bim=bim,fam=fam,pl=pl)
}


## use this instead smarter function - does not use as much RAM
plinkV3<-function(plinkFile){
    pl <- snpStats::read.plink(plinkFile)
    ind<-pl$fam[,1]
    snp<-pl$map[,2]
    geno<-as.integer(as.integer(pl$genotypes)-1)
    dim(geno)<-c(length(ind),length(snp))
    geno[geno==-1]<-NA
    rownames(geno)<-ind
    colnames(geno)<-colnames(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    list(geno=geno,bim=bim,fam=fam,pl=pl)
}



###########################3
writePlink<-function(x,info,file,IID,FID,SEX,PHE){
  n<-ncol(x)
  if(missing(IID))
    IID<-1:n
  if(missing(FID))
    FID<-IID
  if(missing(SEX))
    SEX<-rep(0,n)
  if(missing(PHE))
    PHE<-rep(-9,n)
  
  for(g in GENO2)
    x[x==g]<-paste(strsplit(g,"")[[1]],collapse=" ")

  if(any(is.na(info)))
    stop("no missing in info allowed")
  write.table(cbind(info,x),file=file,col=F,row=F,qu=F,na="0 0")
  write.table(cbind(IID,FID,0,0,SEX,PHE),file=sub("tped","tfam",file),col=F,row=F,qu=F)
  
}

write_tped<-function(pl,out){
  ## function for writing pl plink "object" (read via plinkV2) as .tped and .tfam file - to be able to convert to .bed via plink2  
  print("Please be aware that alleles might get shifted and thus also the genotypes for sites with MAF=0.5")
  print("Also note that one of the alleles might get recoded as '0' for sites with MAF=0.0")
    
  g2<-cbind(pl$bim,t(pl$geno))
  tped<-t(apply(g2,1, function(x) c(paste(x[5],x[5]),paste(x[5],x[6]),paste(x[6],x[6]))[as.numeric(x[7:length(x)])+1]))
  tped[is.na(tped)]<-"0 0"
  g3<-data.frame(g2[,1:4],tped,stringsAsFactors = F)
  
  write.table(g3,paste0(out,".tped"),col=F,qu=F,row=F)
  write.table(pl$fam,paste0(out,".tfam"),col=F,qu=F,row=F)
}


##########################################
sameStrand<-function(A1,B1,A2,B2,wild=0){
    a1<-A1==wild | A1 == A2 | A1==B2
    b1<-B1==wild | B1 == A2 | B1==B2

    a2<-A2==wild | A2 == A1 | A2==B1
    b2<-B2==wild | B2 == A1 | B2==B1
    (a1&b1) | (a2&b2)
}


####
fed<-function(x,n=6){
    if(is.matrix(x) | is.data.frame(x))
        return(x[1:min(nrow(x),n),1:min(ncol(x),n)])
    else
        return(head(x,n))
}
    


ieatIn<-function(a,b){
    ##k<-paste(bim[,1],bim[,4])%in%paste(gl[,1],gl[,2])
    fun<-function(x){
        aa<-a[x,2]
        bb<-b[b[,1]==a[x[1],1],2]
        aa%in%bb
    }
    ta <- tapply(1:nrow(a),factor(a[,1],levels=unique(a[,1])),fun)
    unlist(ta)
}

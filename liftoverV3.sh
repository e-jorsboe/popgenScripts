#!/bin/bash
## make liftover file
Files=$1
chainFile=$2
liftoverProg=/home/emil/software/liftover/liftOver

if [ "$#" -eq 0 ]; then
    echo "supply a plink prefix"
    exit
fi

echo Going to use this chain: $chainFile

if [ -f $Files.bed ]
then
    echo using file $Files.bed
else
    echo the file $Files.bed does not exists  \(use plink --make-bed\)
fi

## make bed file for liftover
Rscript -e "options(scipen=10);map<-read.table('$Files.bim',hea=F,as.is=T);chrs=map[,1];chrs=sub('23','X',chrs);chrs=sub('24','Y',chrs);res<-cbind(paste('chr',chrs,sep=''),map[,4]-1,map[,4],map[,2],0,'+');write.table(res,file='$Files.bedfile',col=F,row=F,qu=F,sep='\t')"

## liftover
$liftoverProg $Files.bedfile $chainFile ${Files}.lifted.bedfile ${Files}.lifted.unmapped

## rm unmapped
grep -v \# ${Files}.lifted.unmapped  | cut -f4 > rmSNPs.txt
## rm duplicates
cat ${Files}.lifted.bedfile | sort -k1 | rev | uniq -d -f3 | rev | cut -f4 > dubSNP.txt
echo removing `wc -l dubSNP.txt`
cat dubSNP.txt >> rmSNPs.txt

## rm duplicates and unmapped
plink --noweb --bfile $Files --exclude rmSNPs.txt --recode --out ${Files}Red 
## change pos

Rscript -e "liftmap<-read.table('${Files}.lifted.bedfile',hea=F,colC=c('character','integer')[c(1,2,2,1,2,1)]);map<-read.table('${Files}Red.map',hea=F,as.is=T); liftmap<-liftmap[liftmap[,4]%in%map[,2],];rownames(map)<-map[,2];map[liftmap[,4],4]<-liftmap[,3];map[liftmap[,4],1]<-sub('chr','',liftmap[,1]);write.table(map,file='hg19.map',col=F,row=F,qu=F)"

echo The new plink file is called  ${Files}.lifted and are on the same strand as the original file
## flip strands 
cat ${Files}.lifted.bedfile | grep -w \- | cut -f4 > flipSNPs.txt
plink  --noweb --ped ${Files}Red.ped --map hg19.map --flip flipSNPs.txt --recode --make-bed --out ${Files}.lifted --allow-extra-chr
## $plink  --noweb --bfile ${Files}Hg19 --recode --out ${Files}Hg19

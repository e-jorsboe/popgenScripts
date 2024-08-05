RF=/home/emil/software/refFinder/refFinder
REF=/emc/emil/viking/data/generated/hg19/merged/hg19NoChr.fa
BIM=$1
echo "-> "
echo "->    looking through reference - lots of ugly non-error messages"
echo "-> "
sed 's/ /\t/g' $BIM | tr -d "\r" | cut -f1,4 |  $RF $REF full | tr '[:lower:]' '[:upper:]' > tmp.ref
echo "-> "
echo "->      done - now the good stuff"
echo "-> "
cat $BIM > tmp.bim


## does work now with [] indexing instead of $
Rscript -e "ref<-read.table('tmp.ref',as.is=T,head=F);bim<-read.table('tmp.bim',as.is=T,head=F);bad<-(ref[,3]!=bim[,5] & ref[,3]!=bim[,6]);write.table(bim[bad,],'bad.bim',col=F,row=F,quote=F)"


echo "#sites Allele1 allele2 ref"

## here sites that have either 0 or N as any of the alleles or in the ref are removed - these might not be put in bad.bim
## as one allele (like 0 A) might agree with ref, or ref might be N
paste $BIM tmp.ref | sed 's/ /\t/g'  | cut -f 5,6,9 | egrep -w -v "0|N" | sort | uniq -c | sort -n -k1 > tmp.res
cat tmp.res
echo "Sites consistent with hg19 plus:"
Rscript -e "r<-read.table('tmp.res',fill=T,as.is=T,head=F);keep<-apply(r[,-1],1,function(x) all(x%in%c('A','C','G','T')));r<-r[keep,];tapply(r[,1],r[,2]==r[,4]|r[,3]==r[,4],sum)"

## so basically tapply loops over all values in r[,1]
## then check if one of bim alleles is in concordance with ref
## enhanced version!!

## diff with anders' because he removes apolymorphic sites - perhaps I kind of should do that too!
## but if aploymorphic and differs from plus - pretty good indication that it is on minus eh?


## make one with bad ones where neither of bim alleles is ref allele

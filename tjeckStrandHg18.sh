RF=/home/emil/software/refFinder/refFinder
REF=/emc/emil/viking/data/generated/hg18_reference/noChr/all.fa
BIM=$1
echo "-> "
echo "->    looking through reference - lots of ugly non-error messages"
echo "-> "
sed 's/ /\t/g' $BIM | tr -d "\r" | cut -f1,4 |  $RF $REF full | tr '[:lower:]' '[:upper:]' > tmp.ref
echo "-> "
echo "->      done - now the good stuff"
echo "-> "
echo "#sites Allele1 allele2 ref"
paste $BIM tmp.ref | sed 's/ /\t/g'  | cut -f 5,6,9 | egrep -w -v "0|N" | sort | uniq -c | sort -n -k1 > tmp.res
cat tmp.res
echo "Sites consistent with hg18 plus:"
Rscript -e "r<-read.table('tmp.res',fill=T,as.is=T,head=F);keep<-apply(r[,-1],1,function(x) all(x%in%c('A','C','G','T')));r<-r[keep,];tapply(r[,1],r[,2]==r[,4]|r[,3]==r[,4],sum)"

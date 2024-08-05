#!/bin/bash

NGSA=/home/emil/git/fastNGSadmix/ngsadmix32

file=$1 ## ONLY filename including extension - be in folder of file when running analysis!
num=$2 ##number of runs per K
P=$3 ##number of threads
out=$4 ##output directory (no /)
K=$5 ##number of populations

##bfile=`basename $file`
mkdir $out

touch $out/$file.$K.likes
rm $out/$file.$K.likes

for f in `seq 1 $num`
do
    echo -n -e $f"\t"
    $NGSA -likes $file -seed $f -K $K -P $P -o $out/$file.$K.$f -printInfo 1
    ##mv $file.$K.$seed.$f.qopt $out/$file.$K.$seed.$f.qopt
    ##mv $file.$K.$seed.$f.fopt $out/$file.$K.$seed.$f.fopt
    ##gzip $out/$file.$K.$seed.$f.fopt.gz
    grep "like=" $out/$file.$K.$f.log | cut -f2 -d " " | cut -f2 -d "=" >> $out/$file.$K.likes
    CONV=`Rscript -e "r<-read.table('$out/$file.$K.likes');r<-r[order(-r[,1]),];cat(sum(r[1]-r<5),'\n')"`

    if [ $CONV -gt 2 ]
    then
        cp $out/$file.$K.$f.qopt $out/$file.$K.$f.qopt_conv
        cp $out/$file.$K.$f.fopt.gz $out/$file.$K.$f.fopt_conv.gz
        cp $out/$file.$K.$f.log $out/$file.$K.$f.log_conv
        break
    fi

done

cat $out/$file.$K.likes | sort -k2 -n -r > $out/$file.$K.likes
## example of run
##./NGSadmixConv.sh filename.beagle.gz 10 4 filename 10

##filename.beagle.gz filename with extensions
##10 runs per seed
##4 Threads
##filename outputname
##10 K (number of pops)

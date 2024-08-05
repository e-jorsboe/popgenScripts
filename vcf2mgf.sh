#!/bin/bash

VCFFILE=$1

VCFFILE2=`echo $VCFFILE | sed 's/\.vcf\.gz//g' | sed 's/\.vcf//g'`

echo $VCFFILE2

bcftools convert --tag GP -g $VCFFILE2 --vcf-ids -o $VCFFILE2 $VCFFILE

ONE=`cat ${VCFFILE2}.samples | wc -l`
TWO=2
INDIS=$(($ONE-$TWO))

pigz -d ${VCFFILE2}.gen.gz 

## SNP-ID, minor-allele, major-allele, dosage in direction of minor-allele (ADDITIVELY)
awk -v s=$INDIS '{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) if( ($(i*3+3)) == "nan" ){ printf ",NA", i; } else{ printf ",%f", $(i*3+3)*2+$(i*3+4) } printf "\n"; }' ${VCFFILE2}.gen > ${VCFFILE2}ADDITIVE.mgf

## SNP-ID, minor-allele, major-allele, dosage in direction of minor-allele (RECESSIVELY)
awk -v s=$INDIS '{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) if( ($(i*3+3)) == "nan" ){ printf ",NA", i; } else{ printf ",%f", $(i*3+3) } printf "\n"; }' ${VCFFILE2}.gen > ${VCFFILE2}RECESSIVE.mgf

## .ann SNP-ID, pos, chr
paste -d" "  <( awk -F " " '{ print $2",", $3"," }' ${VCFFILE2}.gen ) <( cut -f1 -d":" ${VCFFILE2}.gen ) > ${VCFFILE2}.ann

pigz ${VCFFILE2}.gen

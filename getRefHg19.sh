#!/bin/bash

## file should have 2 columns with chr and pos
FILE=$1

RF=/home/emil/software/refFinder/refFinder;
REF=/emc/emil/viking/data/generated/hg19/merged/hg19NoChr.fa;

cat $FILE | tr -d "\r" | $RF $REF full | tr '[:lower:]' '[:upper:]' > ${FILE}.ref

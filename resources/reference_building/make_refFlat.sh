#!/bin/bash

GTF=$1
REF=$2

gtfToGenePred $GTF $REF.tmp
awk '{print $1 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' $REF.tmp > $REF
rm $REF.tmp

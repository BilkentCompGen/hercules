#!/bin/sh
ALIGNMENT=$1 #alignment file
OUTPUT=$2 #output file. note that the output file will be in BAM format
THREAD=$3 #how many thread to use?
MEM=$4 #4G , 800M ...

samtools view -u -b $ALIGNMENT | samtools sort -l 0 -@ $THREAD -m $MEM | samtools rmdup -S - $OUTPUT

#!/bin/sh
LONG=$1
SHORT=$2
OUTFOLDER=$3
THREAD=$4

mkdir -p $OUTFOLDER; bowtie2-build --threads $THREAD -q $LONG $LONG; /usr/bin/time -v -p -o $OUTFOLDER"/alignment.time" bowtie2 --end-to-end -a -f -L 15 --mp 1,1 --np 1 --rdg 0,1 --rfg 0,1 --score-min L,0,-0.12 --no-unal -p $THREAD -x $LONG -U $SHORT | samtools sort -m 8G -l 0 | samtools rmdup -S - $OUTFOLDER"/alignment.bam"


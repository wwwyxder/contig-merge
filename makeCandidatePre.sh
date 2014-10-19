#!/bin/bash
if ! [ -f localTandem.fa ] || ! [ -f localReads.fa ] ; then
    echo 'localTandem.fa or localReads.fa no exist.'
    exit
fi
if ! [ -d bowtie2-2.2.1 ] ; then
    echo "bowtie2 no exist"
    exit
fi
PATH=bowtie2-2.2.1/:$PATH
bowtie2-build localTandem.fa localindex
bowtie2 -f -x localindex -L 30 -N 1 -a -U localReads.fa -S localout.sam --no-head --norc
gawk 'BEGIN {OFS="\t"} $0!~/^@/ {if(and($2,4)==0 ) print $1,$4}' localout.sam |sed 's/:/\t/g' |\
gawk 'BEGIN {OFS="\t"} {if($4+30<$3 || $4>=$2+$3) print}' > local.mid
rm -f localindex.[1234]* localindex.rev.[12]*
echo 'success'

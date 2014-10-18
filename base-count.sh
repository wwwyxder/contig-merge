#!/bin/bash
if [ $# -ne 1 ] ;then
    echo "usage ： $0 <file>"
    exit
fi
cat $1|grep -v "^>"|tr -d "\n"|tr -d "\r"|sed 's/[^ACGT]/[&]/g' > tmp_$1

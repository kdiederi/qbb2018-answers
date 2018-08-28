#!/bin/bash


#to use for other files, need to update ~/qbb2018-answers/day1-homework/SRR072893/SRR072893_map.sam for where new file is
#can change output file names
#numAlign.sam is initial grep (gets rid of useless lines)
#FinalCalc.sam column 1 is count number, then column 2 is what chromosome

grep -v "^@" ~/qbb2018-answers/day1-homework/SRR072893/SRR072893_map.sam | grep -v 2110000 > numAlign.sam

cut -f 3 numAlign.sam | sort | uniq -c > FinalCalc.sam  


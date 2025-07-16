# find the 2 nearest non-overlapping reads for every annotation
bedtools closest -a <(sort -k1,1 -k2,2n $1) -b <(sort -k1,1 -k2,2n $2) -k 2 -io -d -t all 2>/dev/null
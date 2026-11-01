#!/bin/bash
M=$1
cd "$(dirname "$0")"
while read d L S; do
  ./sim_spanning_thresholds $d $L $S $((1000*d+L)) $M > ../data/${M}_d${d}_L${L}.txt
done <<LIST
2 64 4000
2 128 4000
2 256 2000
2 512 1000
2 1024 600
3 16 4000
3 32 2000
3 64 1000
3 128 400
4 8 4000
4 16 2000
4 32 800
6 6 4000
6 12 600
7 4 4000
7 8 600
8 3 4000
8 6 800
5 6 4000
5 12 2000
5 24 250
LIST
echo done > ../data/${M}_finished

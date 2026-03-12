# #!/bin/bash

# 単に値を-1,1に変換
awk 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) $i=($i>0)?1:($i<0)?-1:0; print}' $1 > $2
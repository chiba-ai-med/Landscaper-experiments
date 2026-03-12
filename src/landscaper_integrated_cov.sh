#!/bin/bash
#$ -l nc=4
#$ -p -50
#$ -r yes
#$ -q node.q

#SBATCH -n 4
#SBATCH --nice=50
#SBATCH --requeue
#SBATCH -p node03-06
SLURM_RESTART_COUNT=2

outdir=`echo $5 | sed -e 's|/plot/ratio_group.png||'`

mkdir -p tmp5
cd tmp5
cp /landscaper .
cp /Snakefile .
cp -rf /src .

./landscaper -i ../$1 -o ../$outdir -g ../$2 -d ../$3 -e ../$4 --cores=4 --memgb=10

rm -rf tmp5

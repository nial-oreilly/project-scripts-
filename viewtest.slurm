#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="viewtest"
#SBATCH --output=viewtest.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=32G   # memory per cpu-core
#SBATCH -p MSC

#activate my conda environment
module load Anaconda3
source activate tools

# date written = 13/02/2024 

start_time=$SECONDS

#Data files
FILES=/data4/msc20700669/raw_reads/samtools/view-test/*.bam

for f in $FILES
do
    echo $f
    base=$(basename $f)
    echo $base
    samtools view -bh -X $f $f".bai" -o ${base:0:17}.subset.bam chr2:237,324,003-237,414,328 
done

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

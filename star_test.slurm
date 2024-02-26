#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="star"
#SBATCH --output=star_test.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=60G   # memory per cpu-core

# activate my conda environment
module load Anaconda3
source activate tools

# date written = 12/02/2024 

start_time=$SECONDS
#FILES
FILES=/data4/msc20700669/raw_reads/star_test

for base in 3532
do
    echo $base

# define R1 fastq filename
    fq1=$FILES/${base}_1.fq.gz
# define R2 fastq filename
    fq2=$FILES/${base}_2.fq.gz

    STAR --runThreadN 8 --genomeDir /data4/msc20700669/raw_reads/star/index2 --readFilesCommand gunzip -c --readFilesIn $fq1 $fq2 --outFileNamePrefix $base"_" --outSAMtype BAM SortedByCoordinate
done


elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

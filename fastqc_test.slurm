#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="fastqc"
#SBATCH --output=fastqc_4723.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20G   # memory per cpu-core

# activate my conda environment
module load Anaconda3
source activate tools

# date written = 02/02/2024 

start_time=$SECONDS

cd /data4/msc20700669/raw_reads
echo "current directory..."
pwd
#mkdir outdir
#files="/data4/msc20700669/raw_reads/rna_seq_files.tsv"
#echo "my files"
#cat $files
# run fastqc on all samples
#fastqc --help
#while IFS=$'\t' read -r -a myArray
#do
# echo "reading file..."
# echo "${myArray[0]}"
# var=$(date)
# echo "$var"
# echo "reading reverse read..."
# echo "${myArray[1]}"
fastqc /data/bdigby/Projects/D_O_Connor/raw_reads/4723_R1_001.fastq.gz /data/bdigby/Projects/D_O_Connor/raw_reads/4723_R2_001.fastq.gz -o /data4/msc20700669/outdir_one/
#fastqc "${myArray[0]}" "${myArray[1]}" -o /data4/msc20700669/outdir_one/
#done < $files

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

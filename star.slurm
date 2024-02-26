#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="star"
#SBATCH --output=star.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=60G   # memory per cpu-core
#SBATCH -p MSC

#activate my conda environment
module load Anaconda3
source activate tools

# date written = 08/02/2024 

start_time=$SECONDS

cd /data4/msc20700669/raw_reads/star
echo "current directory..."
pwd
files="/data4/msc20700669/raw_reads/star/trimmed-files.tsv"
echo "my files"
cat $files
# run star on all samples
# star --help
while IFS=$'\t' read -r -a myArray
do
 echo "reading file..."
 echo "${myArray[0]}"
 var=$(date)
 echo "$var"
 echo "reading reverse read..."
 echo "${myArray[1]}"
 STAR --runThreadN 8 --genomeDir /data4/msc20700669/raw_reads/star/index2 --readFilesCommand gunzip -c --readFilesIn "${myArray[0]}" "${myArray[1]}" --outSAMtype BAM SortedByCoordinate
done < $files

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

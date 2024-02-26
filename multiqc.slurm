#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="multiqc_trim"
#SBATCH --output=multiqc_trim.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20G   # memory per cpu-core
#SBATCH -p MSC

#activate my conda environment
module load Anaconda3
source activate tools

# date written = 06/02/2024 

start_time=$SECONDS

cd /data4/msc20700669/raw_reads/multiqc
echo "current directory..."
pwd
files="/data4/msc20700669/raw_reads/trimmed_alpha"
echo "my files"
cat $files

# run multiqc on all samples
#multiqc --help
 multiqc /data4/msc20700669/raw_reads/trimmed_alpha -o /data4/msc20700669/raw_reads/trimmed_alpha
done < $files

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

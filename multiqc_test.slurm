#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="multiqc"
#SBATCH --output=multiqc_4723.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=20G   # memory per cpu-core

# activate my conda environment
module load Anaconda3
source activate tools

# date written = 02/02/2024 

start_time=$SECONDS

cd /data4/msc20700669/raw_reads/multiqc
echo "current directory..."
pwd
multiqc /data4/msc20700669/outdir_one/* -o /data4/msc20700669/outdir_one/
#done < $files

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

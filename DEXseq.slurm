#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="DEXseq2"
#SBATCH --output=Dexseq2.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=64G   # memory per cpu-core
#SBATCH -p highmem
#SBATCH -c 8
#SBATCH --cpus-per-task=8

cd /data4/msc20700669/raw_reads/dexseq/
module load Anaconda3
conda activate tools
Rscript --vanilla DEXSeq.R 

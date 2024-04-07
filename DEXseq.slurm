#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="DEXseqv4"
#SBATCH --output=Dexseqv4.out
#SBATCH --mail-user=n.oreilly18@nuigalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=32G   # memory per cpu-core
#SBATCH -p highmem
#SBATCH -c 8
#SBATCH --cpus-per-task=8

cd /data4/msc20700669/raw_reads/dexseq/Dexseq_v4
module load Anaconda3
conda activate tools
Rscript --vanilla DEXSeq.R 


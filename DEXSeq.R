#!/usr/bin/Rscript

## ----warning=FALSE------------------------------------------------------------
library(GenomicFeatures)
library(DEXSeq)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(XVector)
## ----Preparing the annotation-------------------------------------------------------------------------

download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz", destfile="gencode.v45.annotation.gtf.gz")
txdb = makeTxDbFromGFF("gencode.v45.annotation.gtf.gz")

file.remove("gencode.v45.annotation.gtf.gz")

print("file removed")

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) =
    sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
print("flattenedAnnotation")
## ----Counting reads------------------------------------------------------------

bam_folder <- "/data4/msc20700669/raw_reads/dexseq/bams"
print("bam folder read")
bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)
print("bam files list complete")
bam_files <- c(bam_files)
bamFilesobject = BamFileList( bam_files, asMates = TRUE )
print("bamFilesobject complete")
seqlevelsStyle(flattenedAnnotation) = "UCSC" 
BPPARAM = MulticoreParam(8)
se = summarizeOverlaps(
     flattenedAnnotation, reads = bamFilesobject, singleEnd=FALSE, fragments=TRUE, ignore.strand=TRUE, BPPARAM = BPPARAM)
print("Success !")
## -----------------------------------------------------------------------------
colData(se)$condition = factor(c("tumour", "normal", "tumour", "normal","tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal", "tumour", "normal" ))
colData(se)$patient = factor(c("patient9", "patient9", "patient10", "patient10", "patient2", "patient2", "patient1", "patient1", "patient3", "patient3", "patient4", "patient4", "patient5", "patient5", "patient11", "patient11", "patient6", "patient6", "patient7", "patient7", "patient8", "patient8", "patient12", "patient12"))
dxd = DEXSeqDataSetFromSE( se, design= ~ sample + exon + patient:exon + condition:exon )
print("dxd Success!")


## ----para1,cache=TRUE,results='hide', eval=FALSE------------------------------
#BPPARAM = MultiCoreParam(8)
dxr = DEXSeq(dxd, BPPARAM = BPPARAM)
saveRDS(dxr, "dxr.rds")
print("Finished!")
#dxd = testForDEU( dxd, BPPARAM=BPPARAM)
#dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)




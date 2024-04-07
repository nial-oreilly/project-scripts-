#!/usr/bin/Rscript

## ----warning=FALSE------------------------------------------------------------
library(GenomicFeatures)
library(DEXSeq)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(XVector)

#args = commandArgs(trailingOnly=TRUE)

#if (length(args)==0) {
 # stop("At least one argument must be supplied (summarizedOverlaps object se)", call.=FALSE)
#} else if (length(args)==1) {
  # default output file
 # dxd <- readRDS(args[1])
#}
## ----Preparing the annotation-------------------------------------------------------------------------

download.file("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz", destfile="gencode.v45.annotation.gtf.gz")
txdb = makeTxDbFromGFF("gencode.v45.annotation.gtf.gz")

file.remove("gencode.v45.annotation.gtf.gz")

print("file removed")

flattenedAnnotation = exonicParts( txdb, linked.to.single.gene.only=TRUE )
names(flattenedAnnotation) = sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)
print("flattenedAnnotation")
## ----Counting reads------------------------------------------------------------

bam_folder <- "/data4/msc20700669/raw_reads/dexseq/bams"
print("bam folder read")
bam_files <- list.files(bam_folder, pattern = "\\.bam$", full.names = TRUE)
print("bam files list complete")
bam_files <- c(bam_files)
bamFilesobject = BamFileList( bam_files, asMates = TRUE )
# add this yieldSize = 2e5 if script cancels, between the 2 args
print("bamFilesobject complete")
seqlevelsStyle(flattenedAnnotation) = "UCSC" 
BPPARAM <- BatchtoolsParam(8)
se = summarizeOverlaps(flattenedAnnotation, bamFilesobject, singleEnd=FALSE, fragments=TRUE, ignore.strand=TRUE, BPPARAM = BPPARAM)
print("Success !")
saveRDS(se, "se.rds")
## -----------------------------------------------------------------------------
colData(se)$condition = factor(c("CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN", "CAF", "TAN" ))
colData(se)$patient = factor(c("patient9", "patient9", "patient10", "patient10", "patient2", "patient2", "patient1", "patient1", "patient3", "patient3", "patient4", "patient4", "patient5", "patient5", "patient11", "patient11", "patient6", "patient6", "patient7", "patient7", "patient8", "patient8", "patient12", "patient12"))
formulaFullModel    =  ~ sample + exon + patient:exon + condition:exon
formulaReducedModel =  ~ sample + exon + patient:exon 
dxd = DEXSeqDataSetFromSE( se, design= ~ sample + exon + patient:exon + condition:exon )
print("DEXSeqDataSeqFromSE finished...")
print("writing DEXSeqDataSetFromSE output...")
saveRDS(dxd, file = "dxd2.Rds")




## ----para1,cache=TRUE,results='hide', eval=FALSE------------------------------
print("running estimateSizeFactors...")
dxd = estimateSizeFactors( dxd )
saveRDS(dxd, file = "dxd_estimateSizeFactors.Rds")
print("running estimateDispersions...")
dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM = BPPARAM )
saveRDS(dxd, file = "dxd_estimateDispersions.Rds")
print("running testForDEU...")
dxd = testForDEU( dxd, reducedModel = formulaReducedModel, fullModel = formulaFullModel, BPPARAM = BPPARAM )
saveRDS(dxd, file = "dxd_testForDEU.Rds")
print("running estimateExonFoldChanges...")
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM = BPPARAM)
saveRDS(dxd, file = "dxd_estimateExonFoldChanges.Rds")
print("extracting DEXSeqResults...")
dxr2 = DEXSeqResults( dxd )
saveRDS(dxr2, "dxr_v4.rds")
print("dxr object created. Writing html report")
DEXSeqHTML( dxr2, FDR=0.05, color=c("#FF000080", "#0000FF80") )





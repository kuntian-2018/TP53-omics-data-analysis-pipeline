## This script will generate countmatrices from the BAM files.
## DESeq2 import functions are utilized to generate the count matrices.

library(DESeq2)
## load the transcript annotation file from UCSC.  Make sure to enter the correct genome version
## load the samtools library for R
library(Rsamtools)
##Using the BamFileList function from the Rsamtools package to integrate the BAM files.
## read all bam files in the p53 wild type untreated group

bamFiles4 <- c("/home/ebakker/p53_wt_untreated/M47PT.bam", 
               "/home/ebakker/p53_wt_untreated/M48PT.bam",
               "/home/ebakker/p53_wt_untreated/M52PT.bam",
               "/home/ebakker/p53_wt_untreated/M53PT.bam",
               "/home/ebakker/p53_wt_untreated/M58PT.bam",
               "/home/ebakker/p53_wt_untreated/M608PT.bam",
               "/home/ebakker/p53_wt_untreated/M61PT.bam",
               "/home/ebakker/p53_wt_untreated/M63PT.bam",
               "/home/ebakker/p53_wt_untreated/M67PT.bam",
               "/home/ebakker/p53_wt_untreated/M68PT.bam",
               "/home/ebakker/p53_wt_untreated/M691PT.bam",
               "/home/ebakker/p53_wt_untreated/M697PT.bam",
               "/home/ebakker/p53_wt_untreated/M699PT.bam",
               "/home/ebakker/p53_wt_untreated/M6PT.bam",
               "/home/ebakker/p53_wt_untreated/M700PT.bam",
               "/home/ebakker/p53_wt_untreated/M701PT.bam",
               "/home/ebakker/p53_wt_untreated/M704PT.bam",
               "/home/ebakker/p53_wt_untreated/M70PT.bam",
               "/home/ebakker/p53_wt_untreated/M710PT.bam",
               "/home/ebakker/p53_wt_untreated/M71PT.bam",
               "/home/ebakker/p53_wt_untreated/M73PT.bam",
               "/home/ebakker/p53_wt_untreated/M75PT.bam",
               "/home/ebakker/p53_wt_untreated/M76PT.bam",
               "/home/ebakker/p53_wt_untreated/M80PT.bam",
               "/home/ebakker/p53_wt_untreated/M82PT.bam",
               "/home/ebakker/p53_wt_untreated/M94PT.bam",
               "/home/ebakker/p53_wt_untreated/M98PT.bam")


## read all bam files in the p53 wild type treated group

##bamFiles4 <- c("/home/ebakker/p53_wt_treated/617PT.bam",
##               "/home/ebakker/p53_wt_treated/618PT.bam",
##               "/home/ebakker/p53_wt_treated/634PT.bam",
##               "/home/ebakker/p53_wt_treated/655PT.bam",
##               "/home/ebakker/p53_wt_treated/M101PT.bam",
##               "/home/ebakker/p53_wt_treated/M614PT.bam",
##               "/home/ebakker/p53_wt_treated/M626PT.bam",
##               "/home/ebakker/p53_wt_treated/M632PT.bam",
##               "/home/ebakker/p53_wt_treated/M637PT.bam",
##               "/home/ebakker/p53_wt_treated/M640PT.bam",
##               "/home/ebakker/p53_wt_treated/M645PT.bam",
##               "/home/ebakker/p53_wt_treated/M663PT.bam",
##               "/home/ebakker/p53_wt_treated/M664PT.bam",
##               "/home/ebakker/p53_wt_treated/M666PT.bam",
##               "/home/ebakker/p53_wt_treated/M668PT.bam",
##               "/home/ebakker/p53_wt_treated/M669PT.bam",
##               "/home/ebakker/p53_wt_treated/M684PT.bam",
##               "/home/ebakker/p53_wt_treated/M685PT.bam",
##               "/home/ebakker/p53_wt_treated/M686PT.bam",
##               "/home/ebakker/p53_wt_treated/M687PT.bam",
##               "/home/ebakker/p53_wt_treated/M689PT.bam",
##               "/home/ebakker/p53_wt_treated/M690PT.bam",
##               "/home/ebakker/p53_wt_treated/M694PT.bam",
##               "/home/ebakker/p53_wt_treated/M703PT.bam",
##               "/home/ebakker/p53_wt_treated/M705PT.bam",
##               "/home/ebakker/p53_wt_treated/M714PT.bam",
##               "/home/ebakker/p53_wt_treated/M719PT.bam")


## read all bam files in the p53 mutant untreated group

##bamFiles4 <- c("/home/ebakker/p53_mu_untreated/602PT.bam", 
##              "/home/ebakker/p53_mu_untreated/667PT.bam",
##              "/home/ebakker/p53_mu_untreated/M17PT.bam",
##              "/home/ebakker/p53_mu_untreated/M20PT.bam",
##              "/home/ebakker/p53_mu_untreated/M37PT.bam",
##              "/home/ebakker/p53_mu_untreated/M43PT.bam",
##              "/home/ebakker/p53_mu_untreated/M49PT.bam",
##              "/home/ebakker/p53_mu_untreated/M50PT.bam",
##              "/home/ebakker/p53_mu_untreated/M57PT.bam",
##              "/home/ebakker/p53_mu_untreated/M60PT.bam",
##              "/home/ebakker/p53_mu_untreated/M61PT.bam",
##              "/home/ebakker/p53_mu_untreated/M62PT.bam",
##              "/home/ebakker/p53_mu_untreated/M69PT.bam",
##              "/home/ebakker/p53_mu_untreated/M77PT.bam",
##              "/home/ebakker/p53_mu_untreated/M8PT.bam",
##              "/home/ebakker/p53_mu_untreated/M97PT.bam",
##              "/home/ebakker/p53_mu_untreated/M99PT.bam")

## read all bam files in the p53 mutant treated group

##bamFiles4 <- c("/home/ebakker/p53_mu_treated/649PT.bam")

bamFiles4 <- BamFileList(bamFiles4, yieldSize=1000000)

## The gene models are defined  for counting reads/fragments.

library(GenomicFeatures)
library(GenomicRanges)

## Current work directory is located in the folder storing this R script and the gtf file.
dir <- getwd()
## The gene model is read from an Ensembl GTF file, using Homo_sapiens.GRCh37.o2.gtf.gz file, in order to make a list of exons grouped by gene for counting read.
## None of our sequences are circular using a 0-length character vector.
gtffile <- file.path(dir,"gtf","Homo_sapiens.GRCh37.82.gtf.gz")
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
## Produce a GRangesList of all the exons grouped by gene. Each element of the list a GRanges object of the exons for a gene.
ebg <- exonsBy(txdb, by="gene")
##count the number of reads aligned to each gene

library(GenomicAlignments)
library("BiocParallel")
register(SerialParam())


##This line will create the  SumarizedExperiment object with counts.
##The SumarizedExperiment object will contains a variety of informatin about the experiment.
##The assay of the object contains the matrix of counts. The rowRanges contains the information about the genomic ranges and the colData contains the information about the samples.
se4 <- summarizeOverlaps(features=ebg, reads=bamFiles4,
                         mode="Union",
                         singleEnd=FALSE,
                         ignore.strand=TRUE,
                         fragments=TRUE )

## The count matrix is retrieved from the SumarizedExperiment object.
countdata4 <- assay(se4)
head(countdata4)

## export the count matrix for the p53 wild type untreated group

write.table(countdata4,file="countTableOfp53_wt_untr.txt",sep="\t",quote=FALSE)

## export the count matrix for the p53 wild type treated group

##write.table(countdata4,file="countTableOfp53_wt_tr.txt",sep="\t",quote=FALSE)

## export the count matrix for the p53 mutant untreated group

##write.table(countdata4,file="countTableOfp53_mu_untr.txt",sep="\t",quote=FALSE)

## export the count matrix for the p53 mutant treated group

##write.table(countdata4,file="countTableOfp53_mu_tr.txt",sep="\t",quote=FALSE)
#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if(!require("CNVfilteR", quietly=TRUE)){
    BiocManager::install("CNVfilteR")
}

library(CNVfilteR)

options(scipen = 999)

# Get arguments.
# args <- commandArgs(trailingOnly = TRUE)
# cnv_file <- args[1]
# vcf_file <- args[2]
cnv_file <- "~/data/project/CNVPipe/analysis/res/cnvfilter/sample1.bed"
vcf_file <- "~/data/project/CNVPipe/analysis/res/cnvfilter/sample1.snp.vcf"

# Load copy number data
cnv_gr <- loadCNVcalls(cnvs.file = cnv_file, chr.column = 'chromosome', start.column = 'start', end.column = 'end',
                       cnv.column = 'cnv_type', sample.column = 'sample', genome = 'hg38')

# Load variant data
vcfs <- loadVCFs(vcf.files = vcf_file, cnvs.gr = cnv_gr, min.total.depth = 10)

# Filter
cnv_filter <- filterCNVs(cnv_gr, vcfs)


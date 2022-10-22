#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if(!require("CNVfilteR", quietly=TRUE)){
    BiocManager::install("CNVfilteR")
}

if(!require("dplyr", quietly=TRUE)){
  install.packages('dplyr')
}

library(CNVfilteR)

options(scipen = 999)

# Get arguments.
args <- commandArgs(trailingOnly = TRUE)
cnv_file <- args[1]
vcf_file <- args[2]
out_file <- args[3]

# cnv_file <- "~/data/project/CNVPipe/analysis/res/merge/sample1.bed"
# vcf_file <- "~/data/project/CNVPipe/analysis/snps/freebayes/sample1.snp.vcf"

# Load copy number data
cnv_gr <- loadCNVcalls(cnvs.file = cnv_file, chr.column = 'chromosome', start.column = 'start', end.column = 'end',
                       cnv.column = 'cnv', sample.column = 'sample', genome = 'hg38')
temp_cnv_gr <- trim(cnv_gr)

# Load variant data, only 'PASS' variant will be included
vcfs <- loadVCFs(vcf.files = vcf_file, cnvs.gr = temp_cnv_gr, min.total.depth = 1, genome = 'hg38')

# Filter
cnv_filter <- filterCNVs(temp_cnv_gr, vcfs)

# Get filtered CNVs
filtered <- cnv_filter$cnvs[cnv_filter$cnvs$filter == TRUE]
filtered0 <- data.frame(chromosome = seqnames(filtered),
                        start = start(filtered),
                        end = end(filtered),
                        CNVfilteR = rep('False', length(filtered)))

# Label filtered CNVs in original bed file
cnv <- read.delim(cnv_file)
cnv <- dplyr::left_join(cnv, filtered0, by = c('chromosome', 'start', 'end'))
cnv$CNVfilteR[is.na(cnv$CNVfilteR)] <- 'True'

# write.table(cnv, file = "~/data/project/CNVPipe/analysis/res/cnvfilter/sample1.bed", sep = '\t', quote = F, row.names = F)
write.table(cnv, file = out_file, sep = '\t', quote = F, row.names = F)

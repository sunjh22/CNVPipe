#!/usr/bin/env Rscript

options("repos" = c(CRAN="https://mirror-hk.koddos.net/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   BiocManager::install(version = "3.18")
# }

#if(!require("remotes", quietly=TRUE)){
#  BiocManager::install("remotes")
#}

#Sys.setenv(XML_CONFIG="/usr/bin/xml2-config")

if(!require("CNVfilteR", quietly=TRUE)){
  BiocManager::install("CNVfilteR")
  # BiocManager::install("jpuntomarcos/CNVfilteR", version='1.13.2')
}

if(!require("dplyr", quietly=TRUE)){
  install.packages('dplyr')
}

if(!require("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)){
    BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
}

suppressMessages(library(CNVfilteR))

options(scipen = 999)

# Get arguments.
args <- commandArgs(trailingOnly = TRUE)
cnv_file <- args[1]
vcf_file <- args[2]
out_file <- args[3]
vcf_source <- args[4]

# Load copy number data
cnv_gr <- loadCNVcalls(cnvs.file = cnv_file, chr.column = 'chromosome', start.column = 'start', end.column = 'end',
                       cnv.column = 'cnv', sample.column = 'sample', genome = 'hg38')
temp_cnv_gr <- trim(cnv_gr)

# Load variant data, only 'PASS' variant will be included
vcfs <- loadVCFs(vcf.files = vcf_file, cnvs.gr = temp_cnv_gr, min.total.depth = 5, vcf.source = vcf_source, genome = 'hg38')

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

write.table(subset(cnv, select=-c(sample)), file = out_file, sep = '\t', quote = F, row.names = F)

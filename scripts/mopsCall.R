#!/usr/bin/env Rscript

if (!require("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install(version = "3.18")
}

if (!require("cn.mops", quietly = TRUE)){
    BiocManager::install("cn.mops")
}

if(!require("magrittr", quietly=TRUE)){
    install.packages("magrittr", repos = "http://cran.us.r-project.org")
}

suppressMessages(library(cn.mops))
suppressMessages(library(magrittr))

options(scipen = 999)

# Get arguments.
args <- commandArgs(trailingOnly = TRUE)
x_species <- args[1]
result_dir <- args[2]
bin_size <- as.integer(args[3])
threads <- as.integer(args[4])
bam_files <- tail(args, -4)

# Get arguments from snakemake
# result_dir <- snakemake@params[['resDir']]
# bin_size <- as.integer(snakemake@params[['binSize']])
# threads <- as.integer(snakemake@threads)
# bam_files <- snakemake@input[['bam']]

# Drop MT (fails otherwise) and get read counts in bin_size windows. Windows should contain 50-100 reads each.
human_ref <- paste('chr', c(as.character(seq(22)), "X", "Y"), sep = '')
# Below one is for rice genome
# rice_ref <- c('NC_029256.1','NC_029257.1','NC_029258.1','NC_029259.1','NC_029260.1','NC_029261.1','NC_029262.1','NC_029263.1','NC_029264.1','NC_029265.1','NC_029266.1','NC_029267.1')
# The below one is actually for lettuce genome, do not change the name
rice_ref <- c("CM022518.2","CM022519.2","CM022520.2","CM022521.2","CM022522.2","CM022523.2","CM022524.2","CM022525.2","CM022526.2")

if(x_species=='human'){
    bam_data_ranges <- getReadCountsFromBAM(bam_files, refSeqNames = human_ref, WL = bin_size, parallel = threads)
} else if(x_species=='rice'){
    bam_data_ranges <- getReadCountsFromBAM(bam_files, refSeqNames = rice_ref, WL = bin_size, parallel = threads)
} else{
    bam_data_ranges <- getReadCountsFromBAM(bam_files, WL = bin_size, parallel = threads)
}

# seq_names <- if(x_species=='human') human_ref else rice_ref
# bam_data_ranges <- getReadCountsFromBAM(bam_files, refSeqNames = seq_names, WL = bin_size, parallel = threads)

# Call CNVS and calculate integer copy numbers.
results <- cn.mops(bam_data_ranges, parallel = threads) %>% calcIntegerCopyNumbers()
cnvs <- cnvs(results)

# Function to format cn.mops result GRanges to bed-like data.frame.
granges_to_bed <- function(gr) {
    bed <- data.frame(chrom = as.character(seqnames(gr)), chromStart = start(ranges(gr)), chromEnd = end(ranges(gr)), copyNumber = sub("^CN", "", gr$CN), median = gr$median)
    return(bed)
}

# Function to write bed-like data.frame to bed file.
write_bed <- function(bed, file_name) {
    write.table(bed, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Function to write a particular sample in cn.mops result GRanges to bed file.
write_sample_cnvs <- function(sample_name, cnv_results) {
    sample_cnvs <- cnv_results[cnv_results$sampleName == sample_name]
    bed <- granges_to_bed(sample_cnvs)
    file_name <- sub(".bam$", ".temp.bed", sample_name)
    path <- paste0(result_dir, file_name)
    write_bed(bed, path)
}

# Output CNV bed files for all input samples.
sapply(basename(bam_files), write_sample_cnvs, cnv_results = cnvs)
# save.image(file = paste0(result_dir, "cnmops_wgs.RData"))

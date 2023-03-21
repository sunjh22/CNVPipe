#! /usr/bin/env snakemake

# Include config file, get sample and control, print basic information of CNVPipe run
include: "rules/common.smk"

# Build genome index, BWA index and GATK dictionary for reference genome if they are not existed
include: "rules/pre-processing.smk"

# Reads filtering by Fastp
include: "rules/fastp.smk"

# Map reads with BWA
include: "rules/bwamem.smk"

# Automatically determine bin size
include: "rules/autobin.smk"

# CNV calling by 3 read-depth based methods
include: "rules/cnvkit.smk"
include: "rules/cnvpytor.smk"
include: "rules/cnmops.smk"

# CNV calling by 2 read-pair and split-read based methods
include: "rules/smoove.smk"
include: "rules/delly.smk"

# SNP calling by freebayes or gatk based on read depth
# 4k binSize generally equals to 10X read depth
if config['settings']['gatk-snp'] and config['params']['binSize'] < 4000:
    include: "rules/gatk.smk"
else:
    include: "rules/freebayes.smk"

# Merge CNVs from different tools and assign quality score based on various metrics
include: "rules/merge.smk"

# Get final CNV table
localrules: all

rule all:
    input:
        expand("res/merge/{sample}.bed", sample = config['global']['sample-names']),
        # "cleaned/multiqc-report.html",

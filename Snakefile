#! /usr/bin/env snakemake

# Get sample and control, print basic information of CNVPipe run
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

# CNV calling by 2 read-pair and split-read methods
include: "rules/smoove.smk"
include: "rules/delly.smk"

include: "rules/freebayes.smk"
include: "rules/gatk.smk"

include: "rules/merge.smk"

# Get final CNV table
localrules: all

rule all:
    input:
        expand("res/merge/{sample}.final.bed", sample = config['global']['sample-names']),
        "cleaned/multiqc-report.html",

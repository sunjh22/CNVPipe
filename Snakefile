#! /usr/bin/env snakemake

# Import config file
configfile: "config.yaml"

# Include config file, get sample and control names, print CNVPipe interface
include: "rules/common.smk"

if not config['params']['bam-input']:
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
include: "rules/gatk.smk"
include: "rules/freebayes.smk"

# Merge CNVs from different tools and assign quality score based on various metrics
include: "rules/merge.smk"

# Plot CNVs and export CNV table
include: "rules/report.smk"

# Get final CNV table
if config['params']['species'] != 'human':
    rule all:
        input:
            expand("res/CNVpipe/{sample}.bed", sample = config['global']['sample-names']),
            "cleaned/multiqc-report.html",
else:
    if config['settings']['recurrent']:
        rule all:
            input:
                "res/report/recurrentCNVs",
                expand("res/report/{sample}", sample = config['global']['sample-names']),
                "cleaned/multiqc-report.html",
    else:
        rule all:
            input:
                expand("res/report/{sample}", sample = config['global']['sample-names']),
                "cleaned/multiqc-report.html",

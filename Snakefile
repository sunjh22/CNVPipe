#! /usr/bin/env snakemake

include: "rules/common.smk"

# localrules: all

rule all:
    input:
        # expand("res/merge/{sample}.duphold.vcf", sample = config['global']['sample-names']),
        expand("res/merge/{sample}.duphold.score.bed", sample = config['global']['sample-names']),
        # expand("snps/gatk/{sample}.vqsr.vcf.gz", sample = config['global']['sample-names']),
        expand("snps/freebayes/{sample}.snp.vcf", sample = config['global']['sample-names']),


#SAMPLES = list(filter(lambda f: str(f).startswith('sample'), SAMPLES))

#TOOLS = set(tool for tool in config["TOOLS"])

include: "rules/pre-processing.smk"

include: "rules/cnvkit.smk"
include: "rules/cnvpytor.smk"
# include: "rules/freec.smk"
include: "rules/cnmops.smk"

include: "rules/freebayes.smk"
include: "rules/gatk.smk"

include: "rules/smoove.smk"
include: "rules/delly.smk"

include: "rules/merge.smk"

#include: "rules/cnvfilter.smk"
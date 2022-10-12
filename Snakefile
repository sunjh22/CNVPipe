#! /usr/bin/env snakemake

include: "rules/common.smk"

# rule cnvs:
# 	input:
# 		cnvs_wgs = expand("5k/results/{tool}/{sample}.bed", tool = config["TOOLS"], sample = SAMPLES)
# 	message: "CNV calling complete"


#SAMPLES = list(filter(lambda f: str(f).startswith('sample'), SAMPLES))

#TOOLS = set(tool for tool in config["TOOLS"])

include: "rules/pre-processing.smk"
include: "rules/cnvkit.smk"
include: "rules/cnvpytor.smk"
include: "rules/freec.smk"
include: "rules/cnmops.smk"
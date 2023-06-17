# =================================================================================================
#     Plot CNVs and export CNV table
# =================================================================================================

# TODO 
# 1. For single sampl: a diagram plot showing all CNVs, a dot plot showing the reads distribution of
# top10 (by default) CNVs with the highest overall score, a dot plot showing the SNP BAF in top10
# CNVs.
# 2. For recurrent CNVs, a dot plot showing the reads distribution of 8 samples in the CNV regions;
# a dot plot showing the SNP BAF in the CNV regions

rule plot_single_sample:
    input:
        rules.CNVPipe_prioritize.output,
    output:
        report(
            directory("res/report/{sample}"),
            patterns=["{chrom}_{start}_{end}_{cnv}.png"],
            caption="../report/plotSingleSampleCNVs.rst",
            category="{sample}"
        ),
    params:
        control = (config["global"]["control-sample-names"][0] if 
                    config["global"]["control-sample-names"] else
                    ""),
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/plotCNVs.py"


rule plot_recurrent_cnvs:
    input:
        rules.recurrent_classifyCNV_convert.output,
    output:
        report(
            directory("res/report/recurrentCNVs"),
            patterns=["{chrom}_{start}_{end}_{cnv}.png"],
            caption="../report/plotRecurrentCNVs.rst",
            category="Recurrent CNVs"
        ),
    params:
        control = (config["global"]["control-sample-names"][0] if 
                    config["global"]["control-sample-names"] else ""),
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/plotRecurrentCNVs.py"
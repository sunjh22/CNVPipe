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
        cnv = rules.CNVPipe_prioritize.output,
        bam = "mapped/{sample}.bam",
        snp = (rules.gatk_applyVQSR.output if config['settings']['gatk-snp'] and 
                config['params']['binSize'] < 4000 else rules.freebayes_filter.output.snp),
    output:
        "res/report/{sample}.pdf",
    params:
        absPath = config['params']['absPath'],
        control = config["global"]["control-sample-names"][:3],
    conda:
        "../envs/report.yaml"
    shell:
        "python {params.absPath}/scripts/reportCNVs.py {input.cnv} {input.snp} {input.bam} {output} {params.control}"


rule plot_recurrent_cnvs:
    input:
        rules.recurrent_classifyCNV_convert.output,
    output:
        "res/report/recurrent-CNVs.pdf",
    params:
        absPath = config['params']['absPath']
    conda:
        "../envs/report.yaml"
    shell:
        "python {absPath}/scripts/report.py {input} {output}"
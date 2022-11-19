# =================================================================================================
#     Filter CNV by CNVfilteR
# =================================================================================================

rule cnvfilter_call:
    input:
        bed = "res/merge/{sample}.duphold.score.bed",
        vcf = "snps/freebayes/{sample}.snp.vcf",
    output:
        "res/cnvfilter/{sample}.bed",
    params:
        absPath = config['params']['absPath']
    log:
        "logs/cnvfilter/{sample}.log"
    shell:
        "Rscript {params.absPath}/scripts/cnvFilter.R {input.bed} {input.vcf} {output} > {log} 2>&1"

localrules: all_cnvfilter

rule all_cnvfilter:
    input:
        expand("res/cnvfilter/{sample}.bed", sample = config['global']['sample-names'])
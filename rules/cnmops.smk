# =================================================================================================
#     CNV calling by cn.MOPS
# =================================================================================================

# The results from mops is already in cnv bed format
rule mops_call:
    input:
        expand("mapped/{sample}.bam.bai", sample=config['global']['sample-names']),
        bam = expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
    output:
        bed = expand("res/mops/{sample}.bed", sample=config['global']['sample-names']),
    params:
        resDir = "res/mops/",
        binSize = config['params']['binSize'],
        absPath = config['params']['absPath'],
    threads:
        config['params']['mops']['threads']
    log:
        "logs/mops/call.log"
    benchmark:
        "benchmarks/mops/call.bench"
    conda:
        "../envs/cnmops.yaml"
    shell:
        "Rscript {params.absPath}/scripts/mopsCall.R {params.resDir} {params.binSize} "
        "{threads} {input.bam} > {log} 2>&1"
    # script:
    #     "../scripts/mopsCall.R"

localrules: all_mops

rule all_mops:
    input:
        expand("res/mops/{sample}.bed", sample=config['global']['sample-names']),

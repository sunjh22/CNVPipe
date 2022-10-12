# =================================================================================================
#     CNV calling by cn.MOPS
# =================================================================================================

rule mops_call:
    input:
        bam = expand("mapped/{sample}.bam", sample=config['global']['all-sample-names']),
    output:
        bed = expand("temp/mops/{sample}.txt", sample=config['global']['all-sample-names']),
    params:
        res_dir = "temp/mops/",
        bin_size = config['params']['bin_size'],
    threads:
        config['params']['mops']['threads']
    log:
        "logs/mops/call.log"
    benchmark:
        "benchmarks/mops/call.bench.log"
    conda:
        "envs/cnmops.yaml"
    shell:
        "Rscript scripts/cnmops_wgs.R {params.res_dir} {params.bin_size} {threads} {input.bam} 2> {log}"


rule all_mops:
    input:
        rules.mops_call.output
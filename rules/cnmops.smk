# =================================================================================================
#     CNV calling by cn.MOPS
# =================================================================================================

rule mops_call:
    input:
        bam = expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
        bai = expand("mapped/{sample}.bam.bai", sample=config['global']['sample-names']),
    output:
        bed = expand("res/mops/{sample}.bed", sample=config['global']['sample-names']),
    params:
        res_dir = "res/mops/",
        bin_size = config['params']['bin_size'],
    threads:
        config['params']['mops']['threads']
    log:
        "logs/mops/call.log"
    benchmark:
        "benchmarks/mops/call.bench.log"
    conda:
        "../envs/cnmops.yaml"
    shell:
        "Rscript ../scripts/cnmops_wgs.R {params.res_dir} {params.bin_size} {threads} {input.bam} 2> {log}"

# the results from mops is already in formated cnv bed



rule all_mops:
    input:
        expand("res/mops/{sample}.bed", sample=config['global']['sample-names']),

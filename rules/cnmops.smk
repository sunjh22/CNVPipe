# =================================================================================================
#     CNV calling by cn.MOPS
# =================================================================================================

rule mops_call:
    input:
        bam = expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
    output:
        bed = expand("temp/mops/{sample}.txt", sample=config['global']['sample-names']),
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
        "../scripts/cnmops.R"
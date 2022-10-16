# =================================================================================================
#     CNV calling by Delly
# =================================================================================================

# Delly support defining the window size by the algorithm, but here we still set the window siez
# by ourselves to make it coordinate with other tools
rule delly_call:
    input:
        "mapped/{sample}.bam",
    output:
        cnv = "temp/delly/{sample}.cnv.bcf",
        cov = "temp/delly/{sample}.cov.gz",
    params:
        window = config['params']['bin_size'],
        ref = config['data']['genome'],
        maptrack = config['data']['delly_map'],
    log:
        "logs/delly/{sample}.call.log"
    benchmark:
        "benchmarks/delly/{sample}.benchmark.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "(delly cnv -u -i {params.window} -g {params.ref} "
        "-m {params.maptrack} -c {output.cov} -o {output.cnv} {input}) 2> {log}"

rule delly_classify:
    input:
        rules.delly_call.output.cnv,
    output:
        "temp/delly/{sample}.filtered.bcf",
    log:
        "logs/delly/{sample}.filter.log"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly classify -f germline -o {output} {input} 1> {log}"

rule delly_convert:
    input:
        rules.delly_classify.output,
    output:
        "res/delly/{sample}.bed",
    conda:
        "../envs/delly.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%RDCN]\n' {input} > {output}"

rule all_delly:
    input:
        expand("res/delly/{sample}.bed", sample=config['global']['sample-names'])
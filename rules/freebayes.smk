# =================================================================================================
#     SNP calling by freebayes
# =================================================================================================

rule freebayes_call:
    input:
        ref=config['data']['genome'],
        reads = get_mapped_reads,
    output:
        snp = "snps/freebayes/{sample}.raw.snp.vcf",
    log:
        "logs/freebayes/{sample}.call.log"
    benchmark:
        "benchmarks/freebayes/{sample}.call.bench.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "(freebayes -f {input.ref} {input.reads} | "
        "grep 'TYPE=snp' > {output.snp}) 2> {log}"

rule freebayes_filter:
    input:
        rules.freebayes_call.output.snp,
    output:
        "snps/freebayes/{sample}.snp.vcf",
    log:
        "logs/freebayes/{sample}.filter.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "(bcftools filter -O v -o {output} -s LOWQUAL \ "
        "-e 'QUAL<10 || FMT/DP <5' --set-GTs {input}) 2> {log}"
    
rule all_freebayes:
    input:
        expand("snps/freebayes/{sample}.snp.vcf", sample = config['globla']['sample-names'])

# =================================================================================================
#     SNP calling by freebayes
# =================================================================================================

rule freebayes_call:
    input:
        ref=config['data']['genome'],
        bam = "mapped/{sample}.bam",
        bai = "mapped/{sample}.bam.bai",
    output:
        "snps/freebayes/{sample}.raw.vcf",
    log:
        "logs/freebayes/{sample}.call.log"
    benchmark:
        "benchmarks/freebayes/{sample}.call.bench.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref} {input.bam} > {output} 2> {log}"

rule freebayes_filter:
    input:
        rules.freebayes_call.output,
    output:
        filtered = "snps/freebayes/{sample}.filtered.vcf",
        snp = "snps/freebayes/{sample}.snp.vcf",
    log:
        "logs/freebayes/{sample}.filter.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "(bcftools filter -O v -o {output.filtered} -s LOWQUAL "
        "-e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs . {input};"
        "bcftools view -v snps {output.filtered} > {output.snp}) 2> {log}"
    
rule all_freebayes:
    input:
        expand("snps/freebayes/{sample}.snp.vcf", sample = config['global']['sample-names'])

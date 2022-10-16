# =================================================================================================
#     SNP calling by bcftools
# =================================================================================================

rule bcftools_call:
    input:
        ref=config['data']['genome'],
        reads = get_mapped_reads,
    output:
        bcf = "snps/bcftools/{sample}.bcf",
        vcf = "snps/bcftools/{sample}.raw.vcf.gz",
    log:
        "logs/bcftools/call.log"
    benchmark:
        "benchmarks/bcftools/call.bench.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(samtools mpileup -go {output.bcf} -f {input.ref} {input.reads}; "
        "bcftools call -vmO z -o {output.vcf} {output.bcf}) 2> {log}"

rule bcftools_filter:
    input:
        rules.bcftools_call.output.vcf,
    output:
        filtered = "snps/bcftools/{sample}.filtered.vcf",
        snp = "snps/bcftools/{sample}.snp.vcf",
    log:
        "logs/bcftools/filter.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "(bcftools filter -O v -o {output.filtered} -s LOWQUAL \ "
        "-e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs {input};"
        "bcftools view -v snps {output.filtered} > {output.snp}) 2> {log}"
    
rule all_bcftools:
    input:
        expand("snps/bcftools/{sample}.snp.vcf", sample = config['global']['sample-names'])

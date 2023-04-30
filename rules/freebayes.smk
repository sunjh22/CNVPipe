# =================================================================================================
#     SNP calling by freebayes
# =================================================================================================

# Use freebayes to call SNPs for low-depth data (< 10x)
# In what depth freebayes could provide good results?
rule freebayes_call:
    input:
        "mapped/{sample}.bam.bai",
        ref = config['data']['genome'],
        bam = "mapped/{sample}.bam",
    output:
        "snps/freebayes/{sample}.raw.vcf",
    log:
        "logs/freebayes/{sample}.call.log"
    benchmark:
        "benchmarks/freebayes/{sample}.call.bench"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref} {input.bam} > {output} 2>{log}"

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
        "bcftools view -v snps {output.filtered} > {output.snp}) > {log} 2>&1"
    
localrules: all_freebayes

rule all_freebayes:
    input:
        expand("snps/freebayes/{sample}.snp.vcf", sample = config['global']['sample-names'])

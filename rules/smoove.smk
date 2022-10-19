# =================================================================================================
#     CNV calling by Smoove
# =================================================================================================

# This is a wrapper of Lumpy, which get discordant and split reads automatically, and very fast.
rule smoove_call:
    input:
        "mapped/{sample}.bam"
    output:
        "temp/smoove/{sample}-smoove.genotyped.vcf.gz"
    params:
        outdir = "temp/smoove/",
        exclude = config['data']['smoove_exclude'],
        ref = config['data']['genome'],
    log:
        "logs/smoove/{sample}.call.log"
    benchmark:
        "benchmarks/smoove/{sample}.benchmark.log"
    conda:
        "../envs/smoove.yaml"
    shell:
        "(smoove call --outdir {params.outdir} --exclude {params.exclude} "
        "--name {wildcards.sample} --fasta {params.ref} -p 1 --genotype {input}) 2> {log}"

# Smoove genotype command is a wrapper of svtyper and duphold, the later one will calculate
# the read depth ratio between CNV region and its flanking region, regions with the same
# GC content, and CNV region in the same chromosome.
rule smoove_genotype:
    input:
        vcf = rules.smoove_call.output,
        bam = "mapped/{sample}.bam",
    output:
        "temp/smoove-genotype/{sample}-smoove.genotyped.vcf.gz"
    params:
        outdir = "temp/smoove-genotype/",
        ref = config['data']['genome'],
    log:
        "logs/smoove/{sample}.genotype.log"
    benchmark:
        "benchmarks/smoove/{sample}.genotype.benchmark.log"
    conda:
        "../envs/smoove.yaml"
    shell:
        "smoove genotype -d -x -p 1 --name {wildcards.sample} --outdir {params.outdir} "
        "--fasta {params.ref} --vcf {input.vcf} {input.bam} 2> {log}"

rule smoove_convert:
    input:
        rules.smoove_genotype.output
    output:
        "res/smoove/{sample}.bed"
    log:
        "logs/smoove/{sample}.convert.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "(bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' {input} | "
        "egrep 'DUP|DEL' > {output}) 2> {log}"

rule all_smoove:
    input:
        expand("res/smoove/{sample}.bed", sample=config['global']['sample-names'])
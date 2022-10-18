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

rule smoove_convert:
    input:
        rules.smoove_call.output
    output:
        "res/smoove/{sample}.bed"
    log:
        "logs/smoove/{sample}.convert.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "(bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' {input} | "
        "egrep 'DUP|DEL' > {output}) 2> {log}"

rule all_smoove:
    input:
        expand("res/smoove/{sample}.bed", sample=config['global']['sample-names'])
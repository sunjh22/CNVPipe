# =================================================================================================
#     CNV calling by Smoove
# =================================================================================================

# This is a wrapper of Lumpy, which get discordant and split reads automatically, and very fast.
rule smoove_call:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
    output:
        "temp/smoove/{sample}-smoove.genotyped.vcf.gz",
    params:
        outdir = "temp/smoove/",
        exclude = config['data']['smoove-exclude'],
        ref = config['data']['genome'],
    log:
        "logs/smoove/{sample}.call.log"
    benchmark:
        "benchmarks/smoove/{sample}.bench"
    conda:
        "../envs/smoove.yaml"
    shell:
        "(smoove call --outdir {params.outdir} --exclude {params.exclude} "
        "--name {wildcards.sample} --fasta {params.ref} -p 1 --genotype {input.bam}) > {log} 2>&1"

# Smoove genotype command is a wrapper of svtyper and duphold, the later one will calculatethe read 
# depth ratio between CNV region and its flanking region, and regions with the same GC content, and 
# regions in the same chromosome.
rule smoove_genotype:
    input:
        vcf = rules.smoove_call.output,
        bam = "mapped/{sample}.bam",
    output:
        "temp/smoove-genotype/{sample}-smoove.genotyped.vcf.gz",
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
        "--fasta {params.ref} --vcf {input.vcf} {input.bam} > {log} 2>&1"

# Extract genomic coordinates, transform SVType to copy number (DEL=1, DUP=3), extract QUAL, DHFFC 
# and DHBFC columns
rule smoove_convert:
    input:
        rules.smoove_genotype.output,
    output:
        "res/smoove/{sample}.bed",
    conda:
        "../envs/freebayes.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' "
        "{input} | egrep 'DUP|DEL' | "
        "awk -v OFS='\t' '$4==\"DEL\" && $7<0.7 {{print $1,$2,$3,1,$5\"|\"$7\"|\"$8}} "
        "$4==\"DUP\" && $8>1.3 {{print $1,$2,$3,3,$5\"|\"$7\"|\"$8}}' > {output}"

rule all_smoove:
    input:
        expand("res/smoove/{sample}.bed", sample=config['global']['sample-names'])
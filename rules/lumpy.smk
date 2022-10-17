# =================================================================================================
#     CNV calling by Lumpy
# =================================================================================================

rule lumpy_get_discordants:
    input:
        "mapped/{sample}.bam",
    output:
        "temp/lumpy/{sample}.discordant.bam"
    threads:
        config['params']['samtools']['threads']
    log:
        "logs/lumpy/{sample}.discordant.log"
    conda:
        "../envs/lumpy.yaml"
    shell:
        "(samtools view -b -F 1294 {input} | "
        "samtools sort -@ {threads} -o {output} -) 2> {log}"

rule lumpy_get_splitters:
    input:
        "mapped/{sample}.bam",
    output:
        "temp/lumpy/{sample}.split.bam"
    threads:
        config['params']['samtools']['threads']
    log:
        "logs/lumpy/{sample}.split.log"
    conda:
        "../envs/lumpy.yaml"
    shell:
        "(samtools view -h {input} | "
        "extractSplitReads_BwaMem -i stdin | "
        "samtools sort -@ {threads} -o {output} -) 2> {log}"

rule lumpy_call:
    input:
        bam = "mapped/{sample}.bam",
        discordant = "temp/lumpy/{sample}.discordant.bam",
        split = "temp/lumpy/{sample}.split.bam",
    output:
        "temp/lumpy/{sample}.vcf"
    log:
        "logs/lumpy/{sample}.call.log"
    benchmark:
        "benchmarks/lumpy/{sample}.benmark.log"
    conda:
        "../envs/lumpy.yaml"
    shell:
        "lumpyexpress -B {input.bam} -D {input.discordant} -S {input.split} -o {output} 2> {log}"

rule lumpy_genotype:
    input:
        vcf = rules.lumpy_call.output,
        bam = "mapped/{sample}.bam",
    output:
        genotype = "temp/lumpy/{sample}.gt.vcf",
        json = "temp/lumpy/{sample}.json",
    log:
        "logs/lumpy/{sample}.genotype.log"
    benchmark:
        "benchmarks/lumpy/{sample}.genotype.benchmark.log"
    conda:
        "../envs/lumpy.yaml"
    shell:
        "svtyper -i {input.vcf} -B {input.bam} -l {output.json} > {output.genotype} 2> {log}"

rule lumpy_convert:
    input:
        rules.lumpy_genotype.output.genotype
    output:
        "res/lumpy/{sample}.bed"
    log:
        "logs/lumpy/{sample}.convert.log"
    conda:
        "../envs/freebayes.yaml"
    shell:
        "(bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' {input} | "
        "egrep 'DUP|DEL' > {output}) 2> {log}"

rule all_lumpy:
    input:
        expand("res/lumpy/{sample}.bed", sample=config['global']['sample-names'])
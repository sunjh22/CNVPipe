# =================================================================================================
#     Use duphold to measure read depth of CNV adjacent regions
# =================================================================================================

# Convert CNV bed list into vcf format
rule convert_bed2vcf:
    input:
        bed = "res/merge/{sample}.merged.bed",
        fai = config['data']['genome'] + ".fai",
    output:
        "res/duphold/{sample}.vcf",
    params:
        absPath = config['params']['absPath']
    shell:
        "python {params.absPath}/scripts/bed2vcf.py {input.bed} {input.fai} {output}"

# Calculate DHFFC and DHBFC for CNV region by duphold
rule duphold_score:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
        vcf = rules.convert_bed2vcf.output,
    output:
        "res/duphold/{sample}.duphold.vcf",
    params:
        genome = config['data']['genome'],
    threads: 8
    log:
        "logs/duphold/{sample}.duphold.log"
    benchmark:
        "benchmarks/duphold/{sample}.duphold.bench"
    conda:
        "../envs/smoove.yaml"
    shell:
        "duphold -t {threads} -v {input.vcf} -b {input.bam} -f {params.genome} -o {output}"

# Score CNV region by duphold results
rule score_byDepth:
    input:
        rules.duphold_score.output,
    output:
        bed = "res/duphold/{sample}.duphold.bed",
        scoreBed = "res/duphold/{sample}.duphold.score.bed",
    params:
        absPath = config['params']['absPath']
    conda:
        "../envs/freebayes.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%CN\t%AS\t%DHFFC\t%DHBFC]\t%INFO/TNa\t%INFO/TN\t%INFO/SAMPLE\n' {input} > {output.bed}; "
        "python {params.absPath}/scripts/scoreDuphold.py {output.bed} {output.scoreBed}"


localrules: all_duphold

rule all_duphold:
    input:
        expand("res/duphold/{sample}.duphold.score.bed", sample = config['global']['sample-names'])
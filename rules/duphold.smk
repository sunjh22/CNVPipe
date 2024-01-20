# =================================================================================================
#     Use duphold to measure read depth of CNV adjacent regions
# =================================================================================================

# Convert CNV bed list into vcf format
rule convert_bed2vcf:
    input:
        bed = "res/merge/{sample}.merged.bed",
        fai = config['data']['genome'] + ".fai",
    output:
        "temp/duphold/{sample}.vcf",
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
        "temp/duphold/{sample}.duphold.vcf",
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

rule duphold_extract:
    input:
        rules.duphold_score.output,
    output:
        "temp/duphold/{sample}.duphold.bed",
    conda:
        "../envs/freebayes.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%CN\t%AS\t%DHFFC\t%DHBFC]\t%INFO/TNa\t%INFO/TN\t%INFO/SAMPLE\n' {input} > {output}"

# Score CNV region by duphold results
rule duphold_convert:
    input:
        rules.duphold_extract.output,
    output:
        "res/duphold/{sample}.bed",
    script:
        "../scripts/scoreDuphold.py"


localrules: all_duphold

rule all_duphold:
    input:
        expand("res/duphold/{sample}.bed", sample = config['global']['sample-names'])
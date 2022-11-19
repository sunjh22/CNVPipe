# =================================================================================================
#     Assign scores to CNV
# =================================================================================================

# Merge CNV calls from all tools, at the same time, score the CNVs by calculating the overlap
# fractions with bad genomic region

rule merge_CNVCall:
    input:
        bed = expand(
            "res/{tool}/{sample}.bed",
            tool = ['cnvkit', 'smoove', 'delly', 'mops', 'cnvpytor', 'freec'],
            allow_missing=True
        ),
        low_map = config['data']['smoove-exclude'],
    output:
        "res/merge/{sample}.bed",
    params:
        absPath = config['params']['absPath']
    log:
        "logs/merge/{sample}.merge.log"
    shell:
        "python {params.absPath}/scripts/mergeCNV.py {input.bed} {input.low_map} {output} "
        ">{log} 2>&1"

rule convert_bed2vcf:
    input:
        bed = rules.merge_CNVCall.output,
        fai = config['data']['genome'] + ".fai",
    output:
        "res/merge/{sample}.vcf",
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
        "res/merge/{sample}.duphold.vcf",
    params:
        genome = config['data']['genome'],
    threads: 8
    log:
        "logs/merge/{sample}.duphold.log"
    benchmark:
        "benchmarks/merge/{sample}.duphold.bench"
    conda:
        "../envs/smoove.yaml"
    shell:
        "duphold -t {threads} -v {input.vcf} -b {input.bam} -f {params.genome} -o {output}"

# Score CNV region by duphold results
rule score_byDepth:
    input:
        rules.duphold_score.output.vcf,
    output:
        bed = "res/merge/{sample}.duphold.bed",
        scoreBed = "res/merge/{sample}.duphold.score.bed",
    params:
        absPath = config['params']['absPath']
    conda:
        "../envs/freebayes.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%CN\t%BS\t%DHFFC\t%DHBFC]\t%INFO/TOOL\t%INFO/SAMPLE\n' {input} > {output.bed}; "
        "python {params.absPath}/scripts/scoreDuphold.py {output.bed} {output.scoreBed}"
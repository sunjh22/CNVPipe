# =================================================================================================
#     Single-cell: Count reads, segment, integrate integer copy number
# =================================================================================================

rule counts:
    input:
        bam = "mapped/{sample}.bam",
        bai = "mapped/{sample}.bam.bai",
    output:
        counts = "counts/{sample}.varbin.txt",
    params:
        absPath = config['params']['absPath'],
        boundary = config['data']['bin-boundary'],
        medianCount = "medianCount.txt",
    conda:
        "../envs/pysam.yaml"
    shell:
        "python {params.absPath}/scripts/countReads.py {input.bam} {params.boundary} {output.counts} {params.medianCount}"

rule segment:
    input:
        rules.counts.output.counts,
    output:
        "segment/{sample}.cna.txt",
    params:
        absPath = config['params']['absPath'],
        boundary = config['data']['bin-boundary'],
        ploidyFile = "ploidyEstimate.txt",
        ploidy = config['params']['ploidy'],
    log:
        "logs/segment/{sample}.rscript.log"
    shell:
        "Rscript {params.absPath}/scripts/cbs.plot.r {params.boundary} {input} {output} {params.ploidyFile} {params.ploidy} 1>{log} 2>&1"

rule integrate:
    input:
        expand("segment/{sample}.cna.txt", sample = config['global']['sample-names']),
    output:
        "cnaMatrix.txt",
    params:
        "medianCount.txt",
    script:
        "../scripts/cnaIntegrate.py"
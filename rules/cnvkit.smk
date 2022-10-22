# =================================================================================================
#     CNV calling by CNVKit
# =================================================================================================

# Estimate optimal resolution for specific depth data, will not be triggered automatically
rule cnvkit_autobin:
    input:
        expand("mapped/{sample}.bam.bai", sample=config['global']['sample-names']),
        bam = expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
    output:
        estibin = "logs/cnvkit/estimate.bin",
    params:
        extra = config['params']['cnvkit']['extra'],
        access = config['data']['access'],
        refflat = config['data']['refflat'],
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py autobin {input.bam} {params.extra} -g {params.access} "
        "--annotate {params.refflat} > {output.estibin} 2>&1"

# Call CNVs in samples in batch mode
rule cnvkit_batch:
    input:
        expand("mapped/{sample}.bam.bai", sample=config['global']['sample-names']),
        expand("mapped/{sample}.bam.bai", sample=config['global']['control-sample-names']),
        sample = expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
        control = expand("mapped/{sample}.bam", sample=config['global']['control-sample-names']),
    output:
        reference = "temp/cnvkit/myFlatReference.cnn",
        cns = expand("temp/cnvkit/{sample}.cns", sample=config['global']['sample-names']),
        cnr = expand("temp/cnvkit/{sample}.cnr", sample=config['global']['sample-names']),
    threads:
        config['params']['cnvkit']['threads']
    params:
        ref = config['data']['genome'],
        access = config['data']['access'],
        refflat = config['data']['refflat'],
        outdir = "temp/cnvkit",
        binSize = config['params']['binSize'],
    log:
        "logs/cnvkit/batch.log"
    benchmark:
        "benchmarks/cnvkit/batch.benchmark"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "(cnvkit.py batch {input.sample} -n {input.control} -m wgs -f {params.ref} "
        "--access {params.access} --target-avg-size {params.binSize} -p {threads} "
        "--annotate {params.refflat} --drop-low-coverage --output-reference {output.reference} "
        "-d {params.outdir}) > {log} 2>&1"

# Segment CNVs based on confidence interval
rule cnvkit_segmetric:
    input:
        cns = "temp/cnvkit/{sample}.cns",
        cnr = "temp/cnvkit/{sample}.cnr",
    output:
        "temp/cnvkit/segmetrics/{sample}.cns",
    params:
        "--ci --pi"
    log:
        "logs/cnvkit/segmetric/{sample}.log"        
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py segmetrics -s {input.cns} {input.cnr} {params} -o {output} > {log} 2>&1"

# Call integer copy number
rule cnvkit_call:
    input:
        rules.cnvkit_segmetric.output,
    output:
        "temp/cnvkit/call/{sample}.cns",
    params:
        "-m clonal"
    log:
        "logs/cnvkit/call/{sample}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py call {input} {params} -o {output} > {log} 2>&1"

# Use awk to extract columns: chromosome, start, end, cn, log2, dpeth, probe and weight.
rule cnvkit_convert:
    input:
        rules.cnvkit_call.output,
    output:
        "res/cnvkit/{sample}.bed",
    shell:
        "cat {input} | awk -v OFS='\t' '$8 != 2{{print $1,$2,$3,$8,$5,$9\"|\"$12\"|\"$13}}' "
        "> {output}"

localrules: all_cnvkit

rule all_cnvkit:
    input:
        expand("res/cnvkit/{sample}.bed", sample=config['global']['sample-names']),


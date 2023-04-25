# =================================================================================================
#     CNV calling by CNVKit
# =================================================================================================

# Call CNVs in samples in batch mode
rule cnvkit_batch:
    input:
        get_sample_bai(config['global']['all-sample-names']),
        sample = get_sample_bam(config['global']['sample-names']),
        control = get_sample_bam(config['global']['control-sample-names']),
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
        "echo 'Bin size used for CNV calling is: ' {params.binSize}; "
        "(cnvkit.py batch {input.sample} -n {input.control} -m wgs -f {params.ref} "
        "--access {params.access} --target-avg-size {params.binSize} -p {threads} "
        "--annotate {params.refflat} --drop-low-coverage --output-reference {output.reference} "
        "-d {params.outdir}) > {log} 2>&1"
        
# Calculate confidence interval for log2 ratio and will be used for later call step
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
        "-m clonal --filter ci"
    log:
        "logs/cnvkit/call/{sample}.log"
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "cnvkit.py call {input} {params} -o {output} > {log} 2>&1"

# Merge segments with the same type of CNVs
rule cnvkit_convert:
    input:
        rules.cnvkit_call.output,
    output:
        "res/cnvkit/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
    shell:
        "python {params.absPath}/scripts/cnvkitConvert.py {input} {output}"

localrules: all_cnvkit

rule all_cnvkit:
    input:
        expand("res/cnvkit/{sample}.bed", sample=config['global']['sample-names']),


# =================================================================================================
#     Mapping with BWA MEM
# =================================================================================================

def get_read_group_tags(wildcards):
    res = ["ID:" + wildcards.sample, "SM:" + wildcards.sample, "PL:" + config['params']['bwamem']['platform']]
    return res

def get_bwa_mem_extra(wildcards):
    rg_tags = "\\t".join(get_read_group_tags(wildcards))
    extra = "-R '@RG\\t" + rg_tags + "' " + config["params"]["bwamem"]["extra"]
    return extra

rule map_reads:
    input:
        reads = get_cleaned_reads,
        idx = multiext(config['data']['genome'], ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
        ref = config['data']['genome'],
    output:
        (
            "mapped/{sample}.raw.bam"
            if config['settings']['keep-intermediate']['bwamem']
            else temp("mapped/{sample}.raw.bam")
        )
    params:
        extra = get_bwa_mem_extra,
        sort_extra = config["params"]["samtools"]["sort"],
    threads:
        config['params']['bwamem']['threads'],
    log:
        log = "logs/bwamem/{sample}.log",
    benchmark:
        "benchmarks/bwamem/{sample}.bench",
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "(bwa mem {params.extra} -t {threads} {input.ref} {input.reads} | "
        "samtools sort {params.sort_extra} -@ {threads} -o {output}) >{log} 2>&1"

rule gatk_markDuplicates:
    input:
        rules.map_reads.output,
    output:
        bam = "mapped/{sample}.bam",
        metric = "temp/gatk/{sample}.metric",
    log:
        "logs/gatk/{sample}.markDuplicates.log"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metric} >{log} 2>&1"

rule samtools_index:
    input:
        "mapped/{sample}.bam",
    output:
        "mapped/{sample}.bam.bai",
    threads:
        config['params']['samtools']['threads'],
    log:
        "logs/samtools/{sample}.index.log",
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "samtools index {input} > {log} 2>&1"

localrules: all_bwamem

rule all_bwamem:
    input:
        expand("mapped/{sample}.bam", sample=config['global']['all-sample-names']),
        expand("mapped/{sample}.bam.bai", sample=config['global']['all-sample-names'])

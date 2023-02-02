# =================================================================================================
#     Mapping with BWA MEM
# =================================================================================================

# Build tags for read group used in BWA alignment
def get_bwa_mem_rg(wildcards):
    rg = ["ID:" + wildcards.sample, "SM:" + wildcards.sample, "PL:" + config['params']['bwamem']['platform']]
    tmp_rg_tags = "\\t".join(rg)
    rg_tags = "-R '@RG\\t" + tmp_rg_tags + "' "
    return rg_tags

# Map reads with BWA MEM
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
        read_group = get_bwa_mem_rg,
    threads:
        config['params']['bwamem']['threads']
    log:
        "logs/bwamem/{sample}.log"
    benchmark:
        "benchmarks/bwamem/{sample}.bench"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "(bwa mem -M {params.read_group} -t {threads} {input.ref} {input.reads} | "
        "samtools sort -@ {threads} -o {output}) >{log} 2>&1"

# Mark duplicates with GATK MarkDuplicates
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

# Index bam files
rule samtools_index:
    input:
        "mapped/{sample}.bam",
    output:
        "mapped/{sample}.bam.bai",
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

def get_sample_bam(samples):
    "Quickly access all sample bam file"
    bam = ["mapped/"+sample+".bam" for sample in samples]
    return bam

def get_sample_bai(samples):
    "Quickly access all sample bam index file"
    bai = ["mapped/"+sample+".bam.bai" for sample in samples]
    return bai

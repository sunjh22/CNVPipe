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
        ("mapped/{sample}.raw.bam" if config['settings']['keep-intermediate']['bwamem']
            else temp("mapped/{sample}.raw.bam"))
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
        "samtools sort -@ 10 -o {output}) >{log} 2>&1"

# 1. Mark duplicates with GATK MarkDuplicates
# 2. Recalibrate base quality, prepare for SNP calling, theoretically this step would not affect
# CNV calling. Only bam files produced in this step will be kept, others will be deleted.
# Set different output depends on whether to do BQSR
if config['settings']['bqsr']:
    rule gatk_markDuplicates:
        input:
            rules.map_reads.output,
        output:
            bam = temp("mapped/{sample}.mkdup.bam"),
            metric = "temp/gatk/{sample}.metric",
        threads: 6
        log:
            "logs/gatk/{sample}.markDuplicates.log"
        conda:
            "../envs/pre-processing.yaml"
        shell:
            "gatk MarkDuplicates --java-options \"-Xms10G -Xmx10G -XX:ParallelGCThreads=6\" "
            "-I {input} -O {output.bam} -M {output.metric} >{log} 2>&1"

    rule gatk_BaseRecalibrator:
        input:
            rules.gatk_markDuplicates.output.bam,
        output:
            "mapped/{sample}.recal.table",
        params:
            ref = config['data']['genome'],
            dbsnp = config['data']['gatk-dbsnp'],
            hg38 = config['data']['gatk-hg38'],
            mills = config['data']['gatk-mills'],
        threads: 2
        log:
            "logs/gatk/{sample}.baseRecali.log"
        conda:
            "../envs/pre-processing.yaml"
        shell:
            "gatk --java-options \"-Xms4G -Xmx4G -XX:ParallelGCThreads=2\" BaseRecalibrator "
            "-R {params.ref} -I {input} "
            "-O {output} --known-sites {params.dbsnp} --known-sites {params.hg38} "
            "--known-sites {params.mills} >{log} 2>&1"

    rule gatk_ApplyBQSR:
        input:
            bam = rules.gatk_markDuplicates.output.bam,
            recal = rules.gatk_BaseRecalibrator.output,
        output:
            "mapped/{sample}.bam",
        params:
            ref = config['data']['genome'],
        threads: 2
        log:
            "logs/gatk/{sample}.bqsr.log"
        conda:
            "../envs/pre-processing.yaml"
        shell:
            "gatk --java-options \"-Xms2G -Xmx2G -XX:ParallelGCThreads=2\" ApplyBQSR -R {params.ref} "
            "-I {input.bam} -bqsr {input.recal} -O {output} >{log} 2>&1"
else:
    rule gatk_markDuplicates:
        input:
            rules.map_reads.output,
        output:
            bam = "mapped/{sample}.bam",
        params:
            metric = "temp/gatk/{sample}.metric",
        threads: 6
        log:
            "logs/gatk/{sample}.markDuplicates.log"
        conda:
            "../envs/pre-processing.yaml"
        shell:
            "gatk MarkDuplicates --java-options \"-Xms10G -Xmx10G -XX:ParallelGCThreads=6\" "
            "-I {input} -O {output.bam} -M {params.metric} >{log} 2>&1"

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

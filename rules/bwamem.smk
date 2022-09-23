# =================================================================================================
#     Mapping
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
        reads=get_cleaned_reads,
        idx = multiext(config['data']['genome'], ".amb", ".ann", ".bwt", ".pac", ".sa"),
        # ref=config['data']['genome'],
        # idx=expand(config['data']['genome'] + '.{ext}', ext = ["amb", "ann", "bwt", "pac", "sa", "fai"])
    output:
        (
            "mapped/{sample}.bam"
            if config['settings']['keep-intermediate']['bwamem']
            else temp("mapped/{sample}.bam")
        )
    params:
        extra=get_bwa_mem_extra,
        sorting="samtools",
        sort_order="coordinate",
        sort_extra=config["params"]["samtools"]["sort"],
    threads:
        config['params']['bwamem']['threads'],
    log:
        "logs/bwamem/{sample}.log",
    benchmark:
        "benchmarks/bwamem/{sample}.bench.log",
    wrapper:
        "v1.14.0/bio/bwa/mem"

rule samtools_index:
    input:
        "mapped/{sample}.bam",
    output:
        "mapped/{sample}.bam.bai",
    threads:
        config['params']['samtools']['threads'],
    log:
        "logs/samtools_index/{sample}.log",
    wrapper:
        "v1.14.0/bio/samtools/index"

rule all_bwamem:
    input:
        expand("mapped/{sample}.bam", sample=config['global']['all-sample-names']),
        expand("mapped/{sample}.bam.bai", sample=config['global']['all-sample-names'])

localrules: all_bwamem

def get_mapped_reads():
    return "mapped/{sample}.bam"
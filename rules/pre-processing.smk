# =================================================================================================
#     Pre-processing
# =================================================================================================

# Build index for reference genome
rule samtools_faidx:
    input:
        config["data"]["genome"]
    output:
        config["data"]["genome"] + ".fai"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "samtools faidx {input}"

# Build bwa index for reference genome
if not config['params']['bwamem']['index']:
    rule bwa_index:
        input:
            config["data"]["genome"]
        output:
            multiext(config['data']['genome'], ".amb", ".ann", ".bwt", ".pac", ".sa")
        log:
            "logs/index/bwa_index.log"
        conda:
            "../envs/pre-processing.yaml"
        shell:
            "bwa index {input} > {log} 2>&1"

# Use gatk to build reference dictionary
rule sequence_dictionary:
    input:
        config["data"]["genome"]
    output:
        genome_dict()
    log:
        "logs/index/gatk_dict.log"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output} > {log} 2>&1"

localrules: all_prep

rule all_prep:
    input:
        ref = config["data"]["genome"],
        ref_idx = multiext(config['data']['genome'], ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
        ref_dict = genome_dict(),

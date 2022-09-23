# Helper function to get the name of the genome dictorary file as expected by GATK

def genome_dict():
    return os.path.splitext(config["data"]["genome"])[0] + ".dict"
    
genome=config["data"]["genome"]
genomename=os.path.basename(config["data"]["genome"])
genomedir=os.path.dirname(config["data"]["genome"])

rule samtools_faidx:
    input:
        genome
    output:
        genome + ".fai"
    log:
        "logs/" + genomename + ".samtools_faidx.log"
    wrapper:
        "0.51.3/bio/samtools/faidx"

if not config['params']['bwamem']['index']:
    rule bwa_index:
        input:
            genome
        output:
            genome + ".amb",
            genome + ".ann",
            genome + ".bwt",
            genome + ".pac",
            genome + ".sa"
        log:
            "logs/" + genomename + ".bwa_index.log"
        params:
            prefix=genome,
            algorithm="bwtsw"
        wrapper:
            "0.51.3/bio/bwa/index"

rule sequence_dictionary:
    input:
        genome
    output:
        genome_dict()
    # params:
    #     base= lambda wc: os.path.splitext(genome)[0],
    log:
        "logs/" + genomename + ".sequence_dictionary.log"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output} > {log} 2>&1"


rule all_prep:
    input:
        ref=genome,
        ref_idcs=expand(
            genome + ".{ext}",
            ext=[ "amb", "ann", "bwt", "pac", "sa", "fai" ]
        ),
        ref_dict=genome_dict(),

localrules: all_prep

# Clean up the variables that we used above
del genome
del genomename
del genomedir

include: "fastp.smk"
include: "bwamem.smk"

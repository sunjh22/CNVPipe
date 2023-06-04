# =================================================================================================
#     CNV calling by cn.MOPS
# =================================================================================================

# cn.MOPS requires at least 6 samples to run as it call CNVs based on a mixed possion model.
# cn.MOPS does not use reference genome, thus do not exclude black regions of the genome.
# cn.MOPS requires some least number of reads in one bin, so for low coverage data, 
# the bin size should be reasonablly large to include enough number of reads. 
# It is written in R and is easy to use. The results from mops is already in cnv bed format.
rule mops_call:
    input:
        get_sample_bai(config['global']['sample-names']),
        bam = get_sample_bam(config['global']['sample-names']),
    output:
        bed = expand("res/mops/{sample}.temp.bed", sample=config['global']['sample-names']),
    params:
        species = config['params']['species'],
        resDir = "res/mops/",
        binSize = config['params']['binSize'],
        absPath = config['params']['absPath'],
    threads:
        config['params']['mops']['threads']
    log:
        "logs/mops/call.log"
    benchmark:
        "benchmarks/mops/call.bench"
    conda:
        "../envs/cnmops.yaml"
    shell:
        "Rscript {params.absPath}/scripts/mopsCall.R {params.species} {params.resDir} {params.binSize} "
        "{threads} {input.bam} > {log} 2>&1"

# cn.MOPS may produce results with adjacent CNVs being the same type of CNV but the CN value is 
# slightly different (for example one has cn 3 while the other has cn 4, and they are separated).
# We will try to merge this kind of CNVs and average the copy number.
rule mops_convert:
    input:
        "res/mops/{sample}.temp.bed",
    output:
        "res/mops/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
    script:
        "../scripts/mopsConvert.py"

localrules: all_mops

rule all_mops:
    input:
        expand("res/mops/{sample}.bed", sample=config['global']['sample-names']),

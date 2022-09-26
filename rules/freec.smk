# =================================================================================================
#     CNV calling by Control-FREEC
# =================================================================================================

# extract top 25 lines and first two columns from genome.fai file to be chrLen file
# otherwise all decoy chromosome will be searched and system will throw
# "segmentation fault (core dumped)" error.

rule prepare_chrLen:
    input:
        config['data']['genome'] + '.fai',
    output:
        "temp/freec/chrom.size.txt"
    shell:
        "cut -f 1,2 {input} | sed -n '1,24p' > {output}"

# prepare config file: chrLenFile, window and outputDir is required
rule prepare_configfile:
    input:
        chrLen = rules.prepare_chrLen.output,
    output:
        "temp/freec/configFile.txt",
    params:
        window = config['params']['bin_size'],
        outputDir = "temp/freec/",
        threads = config['params']['freec']['threads'],
        gcprofile = config['data']['GCprofile']
    shell:
        "echo -e '[general]\nchrLenFile={input.chrLen}\nploidy=2\n"
        "breakPointThreshold=1.2\nwindow={params.window}\n"
        "readCountThreshold=10\nmaxThreads={params.threads}\n"
        "GCcontentProfile={params.gcprofile}\n"
        "outputDir={params.outputDir}\n\n[sample]\ninputFormat=BAM\n"
        "mateOrientation=0' > {output}"

# run freec
rule freec_call:
    input:
        config = rules.prepare_configfile.output,
        bam = "mapped/{sample}.bam",
    output:
        cnv = "temp/freec/{sample}.bam_CNVs",
        info = "temp/freec/{sample}.bam_info.txt",
    threads:
        config['params']['freec']['threads']
    log:
        "logs/freec/{sample}.log"
    benchmark:
        "benchmarks/freec/{sample}.bench.log"
    conda:
        "../envs/freec.yaml"
    shell:
        "freec -conf {input.config} -sample {input.bam} 1>{log}"

rule freec_convert:
    input:
        rules.freec_call.output.cnv,
    output:
        "res/freec/{sample}.bed",
    shell:
        "cut -f 1-4 {input} | awk -v OFS='\t' 'BEGIN{{print \"chromosome\tstart\tend\tcn\"}} "
        "{{print $0}}' > {output}"

rule all_freec:
    input:
        expand("res/freec/{sample}.bed", sample=config['global']['sample-names'])


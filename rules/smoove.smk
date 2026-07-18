# =================================================================================================
#     CNV calling by Smoove
# =================================================================================================

# Smoove is a wrapper of Lumpy, duphold and SVtyper, which gets discordant and split reads 
# automatically and runs very fast.
rule smoove_call:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
    output:
        "temp/smoove/{sample}-smoove.genotyped.vcf.gz",
    params:
        outdir = "temp/smoove/",
        exclude = ("--exclude " + config['data']['smoove-exclude'] 
                    if config['data']['smoove-exclude'] else ""),
        ref = config['data']['genome'],
    log:
        "logs/smoove/{sample}.call.log"
    benchmark:
        "benchmarks/smoove/{sample}.bench"
    conda:
        "../envs/smoove.yaml"
    shell:
        "mkdir -p temp/smoove; "
        "touch {output}; "
        "(smoove call --outdir temp/smoove {params.exclude} "
        "--name {wildcards.sample} --fasta {params.ref} -p 1 --genotype {input.bam}) > {log} 2>&1"

# Extract all CNVs (DUP and DEL). Smoove might produce some contradictary calls for noisy or 
# complicated regions, we remove these regions before downstream filtering.
rule smoove_extract:
    input:
        rules.smoove_call.output,
    output:
        "temp/smoove/{sample}.bed",
    conda:
        "../envs/freebayes.yaml"
    shell:
        "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' {input} | "
        "egrep 'DUP|DEL' | awk -v OFS='\t' '$4==\"DEL\" {{print $1,$2,$3,1,$5}} "
        "$4==\"DUP\" {{print $1,$2,$3,3,$5}}' > {output}"
        
#"python {params.absPath}/scripts/smooveFilter.py {output.tmpBed} {output.bed} > {log} 2>&1"
rule smoove_convert:
    input:
        rules.smoove_extract.output,
    output:
        "res/smoove/{sample}.bed",
    script:
        "../scripts/smooveConvert.py"


localrules: all_smoove

rule all_smoove:
    input:
        expand("res/smoove/{sample}.bed", sample=config['global']['sample-names'])

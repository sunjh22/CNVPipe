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

# --------------------------------------------------------------------------------------------------
# Smoove genotype command is a wrapper of svtyper and duphold, the later one will calculatethe read 
# depth ratio between CNV region and its flanking region, and regions with the same GC content, and 
# regions in the same chromosome.
# Change: move the step of duphold genotyping to after merging results
# rule smoove_genotype:
#     input:
#         vcf = rules.smoove_call.output,
#         bam = "mapped/{sample}.bam",
#     output:
#         "temp/smoove-genotype/{sample}-smoove.genotyped.vcf.gz",
#     params:
#         outdir = "temp/smoove-genotype/",
#         ref = config['data']['genome'],
#     log:
#         "logs/smoove/{sample}.genotype.log"
#     benchmark:
#         "benchmarks/smoove/{sample}.genotype.benchmark.log"
#     conda:
#         "../envs/smoove.yaml"
#     shell:
#         "smoove genotype -d -x -p 1 --name {wildcards.sample} --outdir {params.outdir} "
#         "--fasta {params.ref} --vcf {input.vcf} {input.bam} > {log} 2>&1"

# Extract genomic coordinates, transform SVType to copy number (DEL=1, DUP=3), extract QUAL, DHFFC 
# and DHBFC columns
# rule smoove_convert:
#     input:
#         rules.smoove_call.output,
#     output:
#         "res/smoove/{sample}.bed",
#     conda:
#         "../envs/freebayes.yaml"
#     shell:
#         "bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' "
#         "{input} | egrep 'DUP|DEL' | "
#         "awk -v OFS='\t' '$4==\"DEL\" && $7<0.7 {{print $1,$2,$3,1,$5\"|\"$7\"|\"$8}} "
#         "$4==\"DUP\" && $8>1.3 {{print $1,$2,$3,3,$5\"|\"$7\"|\"$8}}' > {output}"

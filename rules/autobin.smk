# =================================================================================================
#     Determine optimal bin size by read depth
# =================================================================================================

# Estimate optimal CNV calling resolution based on the sample with median number of reads
rule autobinBydepth:
    input:
        expand("mapped/{sample}.bam.bai", sample=config['global']['sample-names']),
        bam = expand("mapped/{sample}.bam", sample=config['global']['sample-names']),
    output:
        "logs/autobin/binsize.txt",
    params:
        access = config['data']['access'],
        absPath = config['params']['absPath'],
    conda:
        "../envs/cnvkit.yaml"
    shell:
        "{params.absPath}/scripts/autobin.py {input.bam} -g {params.access} > {output} 2>&1"

localrules: binSizeCal

rule binSizeCal:
    input:
        "logs/autobin/binsize.txt"

if os.path.isfile("logs/autobin/binsize.txt"):
    with open("logs/autobin/binsize.txt", 'r') as f:
        config['params']['binSize'] = int(f.readline().strip())

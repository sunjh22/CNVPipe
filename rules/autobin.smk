# =================================================================================================
#     Determine optimal bin size for CNV calling by read depth
# =================================================================================================

# Estimate optimal CNV calling resolution based on the sample with median number of reads
# This rule cannot be integrated into the whole pipeline, which means if we want to use bin size
# determined by CNVPipe, we need to run CNVPipe twice, the first step calculates the bin
# size and store it into a file, the second step read that file and get bin size.
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

if os.path.isfile("logs/autobin/binsize.txt"):
    with open("logs/autobin/binsize.txt", 'r') as f:
        config['params']['binSize'] = int(f.readline().strip())

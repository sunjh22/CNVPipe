# =================================================================================================
#     CNV calling by cnvpytor
# =================================================================================================

# def get_chroms1():
#     chroms = ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
#     return " ".join(chroms)

# def get_chroms2():
#     chroms = ['NC_029'+str(i)+'.1' for i in range(256,268)]
#     return " ".join(chroms)

def get_chroms(flag):
    if flag:
        return " ".join(['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY'])
    else:
        return " ".join(['NC_029'+str(i)+'.1' for i in range(256,268)])

# CNVpytor is the Python version of CNVnator.
# CNVpytor call CNVs sample by sample, do not require control sample and reference genome, the
# author claimed that pre-defined GC file has been included in the program, which is user-friendly.
# In our test, CNVpytor is IO-consuming.
rule cnvpytor_call:
    input:
        "mapped/{sample}.bam.bai",
        bam = "mapped/{sample}.bam",
    output:
        pytor = "temp/cnvpytor/{sample}.pytor",
        call = "temp/cnvpytor/{sample}.call",
    params:
        binSize = config['params']['binSize'],
        chroms = get_chroms(config['settings']['species']['human'])
    threads:
        config['params']['cnvpytor']['threads']
    log:
        "logs/cnvpytor/{sample}.call.log"
    benchmark:
        "benchmarks/cnvpytor/{sample}.call.bench"
    conda:
        "../envs/cnvpytor.yaml"
    shell:
        "(cnvpytor -root {output.pytor} -j {threads} -chrom {params.chroms} -rd {input.bam}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -his {params.binSize}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -partition {params.binSize}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -call {params.binSize} > {output.call}) > {log} 2>&1"

# "(cnvpytor -root {output.pytor} -j {threads} -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd {input.bam}; \n"

# Filter CNVs based on evalue, q0 (proportion of low-quality reads), pN (proportion of Ns in CNV) 
# and dG (distance from closest large gap)
# Extract columns: chromosome, start, end, cn, log2, evalue, pN, dG.
rule cnvpytor_convert:
    input:
        call = "temp/cnvpytor/{sample}.call",
    output:
        "res/cnvpytor/{sample}.bed",
    script:
        "../scripts/cnvpytorConvert.py"

localrules: all_cnvpytor

rule all_cnvpytor:
    input:
        expand("res/cnvpytor/{sample}.bed", sample = config['global']['sample-names']),

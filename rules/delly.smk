# =================================================================================================
#     CNV calling by Delly
# =================================================================================================

# Call structural variants first
rule delly_call_sv:
    input:
        "mapped/{sample}.bam",
    output:
        "temp/delly/{sample}.sv.bcf"
    params:
        ref = config['data']['genome'],
        exclude = config['data']['smoove_exclude']
    log:
        "logs/delly/{sample}.callsv.log"
    benchmark:
        "benchmarks/delly/{sample}.callsv.benchmark"
    conda:
        "../envs/delly.yaml"
    shell:
        "delly call -g {params.ref} -x {params.exclude} -o {output} {input} 2> {log}"

# Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`,
# here we set window size to 20k.
rule delly_call_cnv:
    input:
        bam = "mapped/{sample}.bam",
        sv = rules.delly_call_sv.output,
    output:
        cnv = "temp/delly/{sample}.cnv.bcf",
        cov = "temp/delly/{sample}.cov.gz",
    params:
        window = config['params']['bin_size'],
        ref = config['data']['genome'],
        maptrack = config['data']['delly_map'],
    log:
        "logs/delly/{sample}.callcnv.log"
    benchmark:
        "benchmarks/delly/{sample}.callcnv.benchmark"
    conda:
        "../envs/delly.yaml"
    shell:
        "(delly cnv -u -i {params.window} -g {params.ref} -l {input.sv} "
        "-m {params.maptrack} -c {output.cov} -o {output.cnv} {input.bam}) 2> {log}"

# Use duphold to genotype delly results
rule delly_genotype:
    input:
        bcf = rules.delly_call_cnv.output.cnv,
        bam = "mapped/{sample}.bam",
    output:
        "temp/delly/{sample}.duphold.vcf"
    params:
        ref = config['data']['genome'],
    log:
        "logs/delly/{sample}.genotype.log"
    benchmark:
        "benchmarks/delly/{sample}.genotype.benchmark.log"
    conda:
        "../envs/smoove.yaml"
    shell:
        "duphold -v {input.bcf} -b {input.bam} -f {params.ref} -o {output}"

# rule delly_classify:
#     input:
#         rules.delly_call.output.cnv,
#     output:
#         "temp/delly/{sample}.filtered.bcf",
#     log:
#         "logs/delly/{sample}.filter.log"
#     conda:
#         "../envs/delly.yaml"
#     shell:
#         "delly classify -f germline -o {output} {input} 1> {log}"

# Transform bcf file to bed and filter low-quality CNVs
rule delly_convert:
    input:
        rules.delly_genotype.output,
    output:
        "res/delly/{sample}.bed",
    conda:
        "../envs/delly.yaml"
    shell:
        "bcftools query -f '%FILTER\t%CHROM\t%POS\t%INFO/END[\t%CN]\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' {input} | "
        "grep 'PASS' | cut -f 2- > {output}"

rule all_delly:
    input:
        expand("res/delly/{sample}.bed", sample=config['global']['sample-names'])
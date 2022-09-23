#print('------------------------------->>>>>>>>>>>')
rule cnvpytor_call:
    input:
        get_mapped_reads(),
    output:
        pytor="temp/cnvpytor/{sample}.pytor",
        call="temp/cnvpytor/{sample}.call",
    params:
        bin_size=config['params']['cnvpytor']['bin_size']
    threads:
        config['params']['cnvpytor']['threads']
    log:
        "logs/cnvpytor/{sample}.call.log"
    benchmark:
        "benchmarks/cnvpytor/{sample}.call.bench.log"
    conda:
        "../envs/cnvpytor.yaml"
    shell:
        "(cnvpytor -root {output.pytor} -j {threads} -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd {input}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -his {params.bin_size}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -partition {params.bin_size}; \n"
        "cnvpytor -root {output.pytor} -j {threads} -call {params.bin_size} > {output.call}) 2> {log}"


# use simple shell commandline to extract columns: chromosome, start, end
# cn, log2, evalue, pN, dG.
rule cnvpytor_convert:
    input:
        call="temp/cnvpytor/{sample}.call",
    output:
        "res/cnvpytor/{sample}.bed",
    script:
        "../scripts/cnvpytor_convert.py"

rule all_cnvpytor:
    input:
        # expand("temp/cnvpytor/{sample}.call", sample=config['global']['sample-names']),
        expand("res/cnvpytor/{sample}.bed", sample=config['global']['sample-names']),

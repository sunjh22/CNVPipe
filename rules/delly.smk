# =================================================================================================
#     CNV calling by Delly
# =================================================================================================

# Delly is written in C++ and calls CNV sample by sample.
# 'delly cnv' optionally take the SV set obtained by 'delly call' as input to improve the accuracy 
# of CNV calling. However 'delly call' is very time-consuming, a 10x sample needs around 24h to 
# complete, so currently we skipped this step and directly used 'delly cnv' for CNV calling.
# Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`.
# Depending on different species, use different Delly calling mode.
if config['params']['species']=='human':
    if config['params']['delly']['complex-mode']:
        rule delly_call_sv_human:
            input:
                "mapped/{sample}.bam.bai",
                bam = "mapped/{sample}.bam",
            output:
                "temp/delly/{sample}.sv.bcf",
            params:
                ref = config['data']['genome'],
                exclude = ("-x " + config['data']['smoove-exclude'] if config['data']['smoove-exclude'] else ""),
            log:
                "logs/delly/{sample}.callsv.log"
            benchmark:
                "benchmarks/delly/{sample}.callsv.bench"
            conda:
                "../envs/delly.yaml"
            shell:
                "delly call -g {params.ref} {params.exclude} -o {output} {input.bam} > {log} 2>&1"

        # Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`.
        rule delly_call_cnv_complex:
            input:
                "mapped/{sample}.bam.bai",
                bam = "mapped/{sample}.bam",
                sv = rules.delly_call_sv_human.output,
            output:
                cnv = "temp/delly/{sample}.cnv.bcf",
                cov = "temp/delly/{sample}.cov.gz",
            params:
                window = config['params']['binSize'],
                ref = config['data']['genome'],
                maptrack = config['data']['delly-map'],
            log:
                "logs/delly/{sample}.callcnv.log"
            benchmark:
                "benchmarks/delly/{sample}.callcnv.benchmark"
            conda:
                "../envs/delly.yaml"
            shell:
                "(delly cnv -u -i {params.window} -g {params.ref} -l {input.sv} "
                "-m {params.maptrack} -c {output.cov} -o {output.cnv} {input.bam}) > {log} 2>&1"
    else:
        rule delly_call_cnv_simple:
            input:
                "mapped/{sample}.bam.bai",
                bam = "mapped/{sample}.bam",
            output:
                cnv = "temp/delly/{sample}.cnv.bcf",
                cov = "temp/delly/{sample}.cov.gz",
            params:
                window = config['params']['binSize'],
                ref = config['data']['genome'],
                maptrack = config['data']['delly-map'],
            log:
                "logs/delly/{sample}.callcnv.log"
            benchmark:
                "benchmarks/delly/{sample}.callcnv.benchmark"
            conda:
                "../envs/delly.yaml"
            shell:
                "(delly cnv -u -i {params.window} -g {params.ref} "
                "-m {params.maptrack} -c {output.cov} -o {output.cnv} {input.bam}) > {log} 2>&1"

    # Extract all CNVs (DUP and DEL) and do the filtering after merging.
    rule delly_extract:
        input:
            "temp/delly/{sample}.cnv.bcf",
        output:
            "temp/delly/extract/{sample}.bed",
        conda:
            "../envs/delly.yaml"
        shell:
            "bcftools query -f '%FILTER\t%CHROM\t%POS\t%INFO/END[\t%CN]\t%QUAL\n' {input} | "
            "grep 'PASS' | cut -f 2- | awk '$4 !=2 {{print$0}}' > {output}"
    
    rule delly_convert:
        input:
            rules.delly_extract.output,
        output:
            "res/delly/{sample}.bed",
        script:
            "../scripts/dellyConvert.py"
else:
    # Use Delly to call structural variants
    rule delly_call_sv_other_species:
        input:
            "mapped/{sample}.bam.bai",
            bam = "mapped/{sample}.bam",
        output:
            "temp/delly/call/{sample}.sv.bcf",
        params:
            ref = config['data']['genome'],
            exclude = ("-x " + config['data']['smoove-exclude'] if config['data']['smoove-exclude'] else ""),
        log:
            "logs/delly/{sample}.callsv.log"
        benchmark:
            "benchmarks/delly/{sample}.callsv.bench"
        conda:
            "../envs/delly.yaml"
        shell:
            "delly call -g {params.ref} {params.exclude} -o {output} {input.bam} > {log} 2>&1"

    # Merge Delly calls
    rule delly_merge:
        input:
            expand("temp/delly/call/{sample}.sv.bcf", sample=config['global']['sample-names']),
        output:
            "temp/delly/mergedSites.sv.bcf",
        log:
            "logs/delly/mergeSV.log"
        conda:
            "../envs/delly.yaml"
        shell:
            "delly merge -o {output} {input} > {log} 2>&1"
    
    # Genotype
    rule delly_genotype:
        input:
            bam = "mapped/{sample}.bam",
            merged = "temp/delly/mergedSites.sv.bcf",
        output:
            "temp/delly/genotype/{sample}.geno.sv.bcf",
        params:
            ref = config['data']['genome'],
            exclude = ("-x " + config['data']['smoove-exclude'] if config['data']['smoove-exclude'] else ""),
        log:
            "logs/delly/{sample}.genotype.log"
        conda:
            "../envs/delly.yaml"
        shell:
            "delly call -g {params.ref} -v {input.merged} -o {output} {params.exclude} {input.bam} > {log} 2>&1"

    # Merge genotyped calls
    rule delly_genotype_merge:
        input:
            expand("temp/delly/genotype/{sample}.geno.sv.bcf", sample=config['global']['sample-names']),
        output:
            "temp/delly/genotype.merged.bcf",
        threads: 10
        log:
            "logs/delly/genotype.merge.log"
        conda:
            "../envs/delly.yaml"
        shell:
            "bcftools merge --threads {threads} -m id -O b -o {output} {input} > {log} 2>&1; "
            "bcftools index {output}"
    
    # Filter calls
    rule delly_filter:
        input:
            "temp/delly/genotype.merged.bcf",
        output:
            "temp/delly/germline.bcf"
        log:
            "logs/delly/filter"
        conda:
            "../envs/delly.yaml"
        shell:
            "delly filter -f germline -o {output} {input} > {log} 2>&1"
    
    # Separate into single samples
    rule delly_uncompress:
        input:
            rules.delly_filter.output,
        output:
            "temp/delly/filter/{sample}.bcf",
        conda:
            "../envs/delly.yaml"
        shell:
            "bcftools view -s {wildcards.sample} -O b -o {output} {input}"

    # Convert bcf file to bed file, in delly result, DEL could have RDCN 2 or even 3,4, DUP could
    # also have RDCN 2, we use strict mode to filter these conflict SVs.
    rule delly_extract:
        input:
            rules.delly_uncompress.output,
        output:
            "temp/delly/extract/{sample}.bed",
        conda:
            "../envs/delly.yaml"
        shell:
            "bcftools query -f '[%FT]\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%RDCN]\n' {input} | "
            "grep 'PASS' | egrep 'DEL|DUP' | cut -f 2- | awk '$4 == \"DEL\" && $5 < 2 {{print$0}} "
            "$4 == \"DUP\" && $5 > 2 {{print$0}}' | cut -f 1,2,3,5 > {output}"
            
    rule delly_convert:
        input:
            rules.delly_extract.output,
        output:
            "res/delly/{sample}.bed",
        script:
            "../scripts/dellyConvert.py"


localrules: all_delly

rule all_delly:
    input:
        expand("res/delly/{sample}.bed", sample=config['global']['sample-names'])

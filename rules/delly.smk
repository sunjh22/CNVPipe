# =================================================================================================
#     CNV calling by Delly
# =================================================================================================

# Delly is written in C++ and calls CNV sample by sample.
# 'delly cnv' optionally take the SV set obtained by 'delly call' as input to improve the accuracy 
# of CNV calling. However 'delly call' is very time-consuming, a 10x sample needs around 24h to 
# complete, so currently we skipped this step and directly used 'delly cnv' for CNV calling.
# Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`.
if config['settings']['species']['human']:
    rule delly_call_cnv:
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
    rule delly_convert:
        input:
            rules.delly_call_cnv.output.cnv,
        output:
            "res/delly/{sample}.bed",
        conda:
            "../envs/delly.yaml"
        shell:
            "bcftools query -f '%FILTER\t%CHROM\t%POS\t%INFO/END[\t%CN]\t%QUAL\n' {input} | "
            "grep 'PASS' | cut -f 2- | awk '$4 !=2 {{print$0}}' > {output}"
else:
    # Use Delly to call structural variants
    rule delly_call_sv:
        input:
            "mapped/{sample}.bam.bai",
            bam = "mapped/{sample}.bam",
        output:
            "temp/delly/call/{sample}.sv.bcf",
        params:
            ref = config['data']['genome'],
            exclude = config['data']['smoove-exclude'],
        log:
            "logs/delly/{sample}.callsv.log"
        benchmark:
            "benchmarks/delly/{sample}.callsv.bench"
        conda:
            "../envs/delly.yaml"
        shell:
            "delly call -g {params.ref} -x {params.exclude} -o {output} {input.bam} > {log} 2>&1"

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
            exclude = config['data']['smoove-exclude'],
        log:
            "logs/delly/{sample}.genotype.log"
        conda:
            "../envs/delly.yaml"
        shell:
            "delly call -g {params.ref} -v {input.merged} -o {output} -x {params.exclude} {input.bam} > {log} 2>&1"

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
    rule delly_convert:
        input:
            rules.delly_uncompress.output,
        output:
            "res/delly/{sample}.bed",
        conda:
            "../envs/delly.yaml"
        shell:
            "bcftools query -f '[%FT]\t%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE[\t%RDCN]\n' {input} | "
            "grep 'PASS' | egrep 'DEL|DUP' | cut -f 2- | awk '$4 == \"DEL\" && $5 < 2 {{print$0}} "
            "$4 == \"DUP\" && $5 > 2 {{print$0}}' | cut -f 1,2,3,5 > {output}"

rule all_delly:
    input:
        expand("res/delly/{sample}.bed", sample=config['global']['sample-names'])

# ---------------------------------------------------------------------------------------------
#! Delly sv call 
# Use Delly to call structural variants first, which could be used as input for CNV calling.
# rule delly_call_sv:
#     input:
#         "mapped/{sample}.bam.bai",
#         bam = "mapped/{sample}.bam",
#     output:
#         "temp/delly/{sample}.sv.bcf",
#     params:
#         ref = config['data']['genome'],
#         exclude = config['data']['smoove-exclude'],
#     log:
#         "logs/delly/{sample}.callsv.log"
#     benchmark:
#         "benchmarks/delly/{sample}.callsv.bench"
#     conda:
#         "../envs/delly.yaml"
#     shell:
#         "delly call -g {params.ref} -x {params.exclude} -o {output} {input.bam} > {log} 2>&1"

# Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`.
# rule delly_call_cnv:
#     input:
#         "mapped/{sample}.bam.bai",
#         bam = "mapped/{sample}.bam",
#         sv = rules.delly_call_sv.output,
#     output:
#         cnv = "temp/delly/{sample}.cnv.bcf",
#         cov = "temp/delly/{sample}.cov.gz",
#     params:
#         window = config['params']['binSize'],
#         ref = config['data']['genome'],
#         maptrack = config['data']['delly-map'],
#     log:
#         "logs/delly/{sample}.callcnv.log"
#     benchmark:
#         "benchmarks/delly/{sample}.callcnv.benchmark"
#     conda:
#         "../envs/delly.yaml"
#     shell:
#         "(delly cnv -u -i {params.window} -g {params.ref} -l {input.sv} "
#         "-m {params.maptrack} -c {output.cov} -o {output.cnv} {input.bam}) > {log} 2>&1"

# Use duphold to genotype delly results
# Change: move the step of duphold genotyping to after merging results
# rule delly_genotype:
#     input:
#         bcf = rules.delly_call_cnv.output.cnv,
#         bam = "mapped/{sample}.bam",
#     output:
#         "temp/delly/{sample}.duphold.vcf",
#     params:
#         ref = config['data']['genome']
#     log:
#         "logs/delly/{sample}.genotype.log"
#     benchmark:
#         "benchmarks/delly/{sample}.genotype.benchmark.log"
#     conda:
#         "../envs/smoove.yaml"
#     shell:
#         "duphold -v {input.bcf} -b {input.bam} -f {params.ref} -o {output} > {log} 2>&1"

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

# Transform bcf file to bed and filter low-quality CNVs by awk instead of delly_classify, because in 
# some cases, value in `FILTER` field is not matched with `FT` tag in `format` field, `FILTER` could 
# be PASS while `FT` could be LowQual, delly_classify will keep this type of CNV, but we want to 
# remove them.
# Extract genomic coordinates, CN and QUAL, DHFFC and DHBFC columns
# rule delly_convert:
#     input:
#         rules.delly_call_cnv.output.cnv,
#     output:
#         "res/delly/{sample}.bed",
#     conda:
#         "../envs/delly.yaml"
#     shell:
#         "bcftools query -f '%FILTER\t%CHROM\t%POS\t%INFO/END[\t%CN]\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' "
#         "{input} | grep 'PASS' | cut -f 2- | "
#         "awk -v OFS='\t' '$4<2 && $7<0.7 {{print $1,$2,$3,$4,$5\"|\"$7\"|\"$8}} "
#         "$4>2 && $8>1.3 {{print $1,$2,$3,$4,$5\"|\"$7\"|\"$8}}'> {output}"

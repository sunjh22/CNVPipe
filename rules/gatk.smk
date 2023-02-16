# =================================================================================================
#     SNP calling by GATK
# =================================================================================================

# We call vcf per sample
# HaplotypeCaller is the most time consuming step
rule gatk_haplotypeCaller:
    input:
        "mapped/{sample}.bam.bai",
        ref_dict = genome_dict(),
        bam = "mapped/{sample}.bam",
    output:
        vcf = "temp/gatk/{sample}.raw_variants.vcf.gz",
    params:
        ref = config['data']['genome'],
        dbsnp = config['data']['gatk-dbsnp'],
        extra = config['params']['gatk']['extra'],
    log:
        "logs/gatk/{sample}.haplotypeCaller.log"
    benchmark:
        "benchmarks/gatk/{sample}.haplotypeCaller.bench"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "gatk {params.extra} HaplotypeCaller -R {params.ref} -I {input.bam} -O {output.vcf} "
        "--dbsnp {params.dbsnp} >{log} 2>&1"

rule gatk_variantRecalibrator:
    input:
        rules.gatk_haplotypeCaller.output.vcf,
    output:
        recal = "temp/gatk/{sample}.recal",
        tranches = "temp/gatk/{sample}.tranches",
        rscript = "temp/gatk/{sample}.plots.R",
    params:
        ref = config['data']['genome'],
        #dbsnp = config['data']['gatk-dbsnp'],
        hapmap = config['data']['gatk-hapmap'],
        omni = config['data']['gatk-omni'],
        geno1000 = config['data']['gatk-1000g'],
        extra = config['params']['gatk']['extra'],
    log:
        "logs/gatk/{sample}.variantRecalibrator.log"
    benchmark:
        "benchmarks/gatk/{sample}.variantRecalibrator.bench"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "gatk {params.extra} VariantRecalibrator -R {params.ref} -V {input} "
        "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} "
        "--resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} "
        "--resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.geno1000} "
        "-an MQ -an MQRankSum -an QD -an ReadPosRankSum -an FS -an SOR -mode SNP "
        "--max-gaussians 4 "
        "-O {output.recal} --tranches-file {output.tranches} --rscript-file {output.rscript} "
        ">{log} 2>&1"

# the variance of these two is zero, may be due to simulation data, should be added for real data.
# -an MQ -an MQRankSum

rule gatk_applyVQSR:
    input:
        vcf = rules.gatk_haplotypeCaller.output.vcf,
        recal = rules.gatk_variantRecalibrator.output.recal,
        tranches = rules.gatk_variantRecalibrator.output.tranches,
    output:
        "snps/gatk/{sample}.vqsr.vcf.gz",
    params:
        ref = config['data']['genome'],
        extra = config['params']['gatk']['extra'],
    log:
        "logs/gatk/{sample}.vqsr.log"
    benchmark:
        "benchmarks/gatk/{sample}.vqsr.bench"
    conda:
        "../envs/pre-processing.yaml"
    shell:
        "gatk {params.extra} ApplyVQSR -R {params.ref} -V {input.vcf} -O {output} "
        "--truth-sensitivity-filter-level 99.9 --create-output-variant-index true "
        "--tranches-file {input.tranches} --recal-file {input.recal} -mode SNP >{log} 2>&1"

localrules: all_gatk

rule all_gatk:
    input:
        expand("snps/gatk/{sample}.vqsr.vcf.gz", sample = config['global']['sample-names']),

# As we just kept SNP in vcf, this selecting step is not necessary.
# rule gatk_selectVariants:
#     input:
#         rules.gatk_applyVQSR.output,
#     output:
#         "snps/gatk/{sample}.vqsr.snp.vcf.gz",
#     params:
#         ref = config['data']['genome'],
#         extra = config['params']['gatk']['extra'],
#     log:
#         "logs/gatk/{sample}.selectVariants.log"
#     benchmark:
#         "benchmarks/gatk/{sample}.selectVariants.bench"
#     conda:
#         "../envs/pre-processing.yaml"
#     shell:
#         "gatk {params.extra} SelectVariants -R {params.ref} -V {input} -O {output} "
#         "--select-type-to-include SNP --exclude-filtered "
#         "--exclude-non-variants --remove-unused-alternates >{log} 2>&1"

# =================================================================================================
#     Assign various scores to CNV
# =================================================================================================

# Merge CNV calls from 5 tools according to the read depth.
# if RD < 1x: merge cn.mops, cnvkit and cnvpytor
# if 1x < RD < 5x: merge cn.mops, cnvkit, delly, cnvpytor and smoove
# if RD > 5x: merge smoove, delly, cnvkit, cnvpytor, mops
# Assign 'accumulative score' (1. AS)
if config['params']['binSize'] < 8000:
    rule merge_CNVCall_highDepth:
        input:
            bed = expand(
                "res/{tool}/{sample}.bed", tool = ['smoove', 'delly', 'cnvkit', 'cnvpytor', 'mops'],
                allow_missing=True
            ),
        output:
            "res/merge/{sample}.merged.bed",
        params:
            absPath = config['params']['absPath']
        log:
            "logs/merge/{sample}.merge.log"
        shell:
            "python {params.absPath}/scripts/mergeCNV.py {input.bed} {output} >{log} 2>&1"
elif 8000 < config['params']['binSize'] < 40000:
    rule merge_CNVCall_medianDepth:
        input:
            bed = expand(
                "res/{tool}/{sample}.bed", tool = ['mops', 'cnvkit', 'delly', 'cnvpytor', 'smoove'],
                allow_missing=True
            ),
        output:
            "res/merge/{sample}.merged.bed",
        params:
            absPath = config['params']['absPath']
        log:
            "logs/merge/{sample}.merge.log"
        shell:
            "python {params.absPath}/scripts/mergeCNVMedianDepth.py {input.bed} {output} >{log} 2>&1"
else:
    rule merge_CNVCall_lowDepth:
        input:
            bed = expand(
                "res/{tool}/{sample}.bed", tool = ['mops', 'cnvkit', 'cnvpytor'],
                allow_missing=True
            ),
        output:
            "res/merge/{sample}.merged.bed",
        params:
            absPath = config['params']['absPath']
        log:
            "logs/merge/{sample}.merge.log"
        shell:
            "python {params.absPath}/scripts/mergeCNVLowDepth.py {input.bed} {output} >{log} 2>&1"

localrules: all_merge_CNVCall
rule all_merge_CNVCall:
    input:
        expand("res/merge/{sample}.merged.bed", sample = config['global']['sample-names'])

# Apply duphold and assign 'depth score' (2. DS)
include: "duphold.smk"

# Apply cnvfilter and assign 'SNP score' (3. SS) if SNPs could be called
if config['settings']['gatk-snp'] and config['params']['binSize'] < 4000:
    rule cnvfilter_call_gatk:
        input:
            bed = rules.score_byDepth.output.scoreBed,
            vcf = rules.gatk_applyVQSR.output,
        output:
            "res/cnvfilter/{sample}.bed",
        params:
            absPath = config['params']['absPath'],
            vcf_source = "HaplotypeCaller",
        log:
            "logs/cnvfilter/{sample}.log"
        conda:
            "../envs/cnvfilter.yaml"
        shell:
            "Rscript {params.absPath}/scripts/cnvFilter.R {input.bed} {input.vcf} {output} {params.vcf_source} > {log} 2>&1"
else:
    rule cnvfilter_call_freebayes:
        input:
            bed = rules.score_byDepth.output.scoreBed,
            vcf = rules.freebayes_filter.output.snp,
        output:
            "res/cnvfilter/{sample}.bed",
        params:
            absPath = config['params']['absPath'],
            vcf_source = "freeBayes",
        log:
            "logs/cnvfilter/{sample}.log"
        conda:
            "../envs/cnvfilter.yaml"
        shell:
            "Rscript {params.absPath}/scripts/cnvFilter.R {input.bed} {input.vcf} {output} {params.vcf_source} > {log} 2>&1"

# Calculate overlap fraction with low-complexity region and assign 'good score' (4. GS)
# Calculate overlap fraction with CNVs in normal population and assign 'normal score' (5. NS)
rule good_normal_score:
    input:
        "res/cnvfilter/{sample}.bed",
    output:
        "res/merge/{sample}.goodscore.bed",
    params:
        absPath = config['params']['absPath'],
        badList = config['data']['smoove-exclude'],
        normalList = config['data']['normal-common-cnv'],
    log:
        "logs/merge/{sample}.goodscore.log"
    shell:
        "python {params.absPath}/scripts/goodNormalScore.py {input} {params.badList} "
        "{params.normalList} {output} >{log} 2>&1"

# Apply ClassifyCNV and assign 'pathogenicity score' (6. PS)
rule classifycnv_predict:
    input:
        rules.good_normal_score.output,
    output:
        "res/classifycnv/{sample}.classifycnv.txt",
    params:
        absPath = config['params']['absPath']
    threads: 2
    log:
        "logs/merge/{sample}.patho.log"
    conda:
        "../envs/classifycnv.yaml"
    shell:
        "python {params.absPath}/scripts/classifyCNV.py --absPath {params.absPath} --infile {input}"
        " --GenomeBuild hg38 --cores {threads} >{log} 2>&1"

rule classifycnv_convert:
    input:
        normal_bed = rules.good_normal_score.output,
        patho_bed = rules.classifycnv_predict.output,
    output:
        "res/merge/{sample}.bed"
    params:
        absPath = config['params']['absPath']
    shell:
        "python {params.absPath}/scripts/classifyCNVConvert.py {input.normal_bed} {input.patho_bed}"
        " {output}"
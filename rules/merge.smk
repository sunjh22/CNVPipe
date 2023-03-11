# =================================================================================================
#     Assign various scores to CNV
# =================================================================================================

# Merge CNV calls from 5 tools, if two CNVs have overlaps, we extend the breakpoints.
# Assign 'accumulative score' (1. AS)
# If binSize is larger than 80k, which means read depth is lower than 0.5x, Delly and Smoove
# will not work, thus we will only merge the results from cnvkit, cnvpytor and cn.mops
if config['params']['binSize'] < 80000:
    rule merge_CNVCall:
        input:
            bed = expand(
                "res/{tool}/{sample}.bed", tool = ['cnvkit', 'delly', 'mops', 'cnvpytor', 'smoove'],
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
else:
    rule merge_CNVCall_lowDepth:
        input:
            bed = expand(
                "res/{tool}/{sample}.bed", tool = ['cnvkit', 'mops', 'cnvpytor'],
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
rule cnvfilter_call:
    input:
        bed = rules.score_byDepth.output.scoreBed,
        vcf = rules.freebayes_filter.output.snp,
    output:
        "res/cnvfilter/{sample}.bed",
    params:
        absPath = config['params']['absPath']
    log:
        "logs/cnvfilter/{sample}.log"
    shell:
        "Rscript {params.absPath}/scripts/cnvFilter.R {input.bed} {input.vcf} {output} > {log} 2>&1"

# Calculate overlap fraction with low-complexity region and assign 'good score' (4. GS)
# Calculate overlap fraction with CNVs in normal population and assign 'normal score' (5. NS)
rule good_normal_score:
    input:
        rules.cnvfilter_call.output,
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
        "res/merge/{sample}.final.bed"
    params:
        absPath = config['params']['absPath']
    shell:
        "python {params.absPath}/scripts/classifyCNVConvert.py {input.normal_bed} {input.patho_bed}"
        " {output}"
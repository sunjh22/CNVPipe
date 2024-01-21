# =================================================================================================
#     Assign various scores to CNV
# =================================================================================================

# Merge CNV calls from 5 tools by different strategies according to the read depth of the data.
# if RD < 1x (bin > 40k): merge cn.mops, cnvkit and cnvpytor
# if 1x <= RD < 5x (8k < bin < 40k): merge cn.mops, cnvkit, delly, cnvpytor and smoove
# if RD >= 5x (bin < 8k): merge smoove, delly, cnvkit, cnvpytor, mops
# Assign 'accumulative score' (1. AS)
if config['params']['binSize'] <= 8000:
    rule merge_CNVCall_highDepth:
        input:
            bed = expand(
                "res/{tool}/{sample}.bed", tool = ['smoove', 'delly', 'cnvkit', 'cnvpytor', 'mops'],
                allow_missing=True
            ),
        output:
            "res/merge/{sample}.merged.bed",
        params:
            flag = 'high'
        script:
            "../scripts/mergeCNV.py"
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
            flag = 'median'
        script:
            "../scripts/mergeCNV.py"
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
            flag = 'low'
        script:
            "../scripts/mergeCNV.py"

localrules: all_merge_CNVCall
rule all_merge_CNVCall:
    input:
        expand("res/merge/{sample}.merged.bed", sample = config['global']['sample-names'])

# Apply duphold and assign 'depth score' (2. DS)
include: "duphold.smk"

# Output CNVpipe results at this step for other species (like rice), following steps will not be ran.
# Note the difference of 'CNVpipe' here and 'CNVPipe' in later rules.
# Selection criteria for CNV is duphold_score>0 and tool_num>=2.
rule cnvpipe_convert:
    input:
        rules.duphold_convert.output,
    output:
        "res/CNVpipe/{sample}.bed",
    shell:
        "awk '$7>0 && $12>=2' {input} | cut -f 1-13 > {output}"

localrules: all_cnvpipe_convert
rule all_cnvpipe_convert:
    input:
        expand("res/CNVpipe/{sample}.bed", sample = config['global']['sample-names'])

# Apply cnvfilter, and assign 'True' for bad CNVs that cannot pass SNP verification, otherwise 
# label 'True'.
rule cnvfilter_call:
    input:
        bed = rules.duphold_convert.output,
        vcf = (rules.gatk_applyVQSR.output if config['settings']['gatk-snp'] and config['params']['binSize'] < 4000
                else rules.freebayes_filter.output.snp),
    output:
        "res/cnvfilter/{sample}.bed",
    params:
        absPath = config['params']['absPath'],
        vcf_source = ("HaplotypeCaller" if config['settings']['gatk-snp'] and config['params']['binSize'] < 4000
                        else "freeBayes"),
    log:
        "logs/cnvfilter/{sample}.log"
    conda:
        "../envs/cnvfilter.yaml"
    shell:
        "Rscript {params.absPath}/scripts/cnvFilter.R {input.bed} {input.vcf} {output} {params.vcf_source} > {log} 2>&1"           

# Calculate overlap fraction with low-complexity region and assign 'good score' (4. GS)
# Calculate overlap fraction with low-mappable regions and assign 'map score' (. MS)
# Calculate overlap fraction with CNVs in normal population and assign 'normal score' (5. NS)
rule good_normal_score:
    input:
        "res/cnvfilter/{sample}.bed",
    output:
        "res/merge/{sample}.goodscore.bed",
    params:
        badList = config['data']['smoove-exclude'],
        lowMapList = config['data']['low-mappable'],
        normalList = config['data']['normal-common-cnv'],
    script:
        "../scripts/goodNormalScore.py"

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

# Convert ClassifyCNV output to well-formated CNV bed file.
rule classifycnv_convert:
    input:
        normal_bed = rules.good_normal_score.output,
        patho_bed = rules.classifycnv_predict.output,
    output:
        "res/classifycnv/{sample}.bed",
    params:
        "before-recurrent"
    script:
        "../scripts/classifyCNVConvert.py"

# Filter CNVs based on toolNum, DS and cnvfilter label.
rule CNVPipe_convert:
    input:
        rules.classifycnv_convert.output,
    output:
        "res/CNVPipe/{sample}.bed",
    shell:
        "(head -n1 {input}; awk '$12>=2 && $7>0 && $14==\"True\"' {input}) > {output}"

# Prioritize CNVs based on all score metrics, pathogenicity has very big weight.
rule CNVPipe_prioritize:
    input:
        rules.CNVPipe_convert.output,
    output:
        "res/CNVPipe/{sample}.priority.bed",
    script:
        "../scripts/cnvpipeConvert.py"

# Optionally identify recurrent CNVs for samples with specific phenotypes
rule find_recurrent_cnvs:
    input:
        expand("res/CNVPipe/{sample}.bed", sample=config['global']['sample-names']),
    output:
        "temp/recurrent/recurrent.bed",
    params:
        sampleNumThe = config['params']['recurrent-threshold'],
    script:
        "../scripts/findRecurrentCNV.py"

rule recurrent_sort:
    input:
        rules.find_recurrent_cnvs.output,
    output:
        "temp/recurrent/recurrent.sort.bed",
    shell:
        "sort -Vk 1 -k 2,3n {input} > {output}"

rule recurrent_classifyCNV:
    input:
        rules.recurrent_sort.output,
    output:
        "res/classifycnv/recurrent.classifycnv.txt",
    params:
        absPath = config['params']['absPath'],
    log:
        "logs/recurrent/recurrent.classifyCNV.log"
    shell:
        "python {params.absPath}/scripts/classifyCNV.py --absPath {params.absPath} --infile "
        "{input} --GenomeBuild hg38 --cores {threads} >{log} 2>&1; "

rule recurrent_classifyCNV_convert:
    input:
        bed = rules.recurrent_sort.output,
        classify = rules.recurrent_classifyCNV.output,
    output:
        "res/recurrent/recurrent.bed",
    params:
        "after-recurrent"
    script:
        "../scripts/classifyCNVConvert.py"
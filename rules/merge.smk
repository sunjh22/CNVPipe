# =================================================================================================
#     Assign various scores to CNV
# =================================================================================================

# Merge CNV calls from 4 tools, if two CNVs have overlaps, we extend the breakpoints.
# Assign 'accumulative score' (1. AS)
rule merge_CNVCall:
    input:
        bed = expand(
            "res/{tool}/{sample}.bed", tool = ['cnvkit', 'delly', 'mops', 'cnvpytor'],
            allow_missing=True
        ),
    output:
        "res/merge/{sample}.bed",
    params:
        absPath = config['params']['absPath']
    log:
        "logs/merge/{sample}.merge.log"
    shell:
        "python {params.absPath}/scripts/mergeCNV3.py {input.bed} {output} >{log} 2>&1"

# Apply duphold and assign 'depth score' (2. DS)
include: "duphold.smk"

# Apply cnvfilter and assign 'SNP score' (3. SS) if SNPs could be called
include: "cnvfilter.smk"

# Calculate overlap fraction with low-complexity region and assign 'good score' (4. GS)
# Calculate overlap fraction with CNVs in normal population and assign 'normal score' (5. NS)
rule good_normal_score:
    input:
        rules.cnvfilter_call.output
    output:
        "res/merge/{sample}.goodscore.bed"
    params:
        absPath = config['params']['absPath']
    shell:
        "python {params.absPath}/scripts/goodNormalScore.py {input.bed} {output} >{log} 2>&1"

# Apply ClassifyCNV and assign 'pathogenicity score' (6. PS)
include: "classifycnv.smk"

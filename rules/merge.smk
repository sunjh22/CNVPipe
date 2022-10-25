# =================================================================================================
#     Merge CNV calling results
# =================================================================================================

rule merge_call:
    input:
        bed = expand(
            "res/{tool}/{sample}.bed",
            tool = ['cnvkit', 'cnvpytor', 'freec', 'mops', 'smoove', 'delly'],
            allow_missing=True
        ),
        low_map = config['data']['smoove-exclude'],
    output:
        "res/merge/{sample}.bed",
    params:
        absPath = config['params']['absPath']
    log:
        "logs/merge/{sample}.merge.log"
    shell:
        "python {params.absPath}/scripts/mergeReadDepthCNV.py {input.bed} {input.low_map} {output} "
        ">{log} 2>&1"

localrules: all_merge

rule all_merge:
    input:
        expand("res/merge/{sample}.bed", sample = config['global']['sample-names'])

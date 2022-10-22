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
    shell:
        "python ../scripts/mergeReadDepthCNV.py {input.bed} {input.low_map} {output}"

localrules: all_merge

rule all_merge:
    input:
        expand("res/merge/{sample}.bed", sample = config['global']['sample-names'])

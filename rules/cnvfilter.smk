# =================================================================================================
#     Score CNV by CNVfilteR
# =================================================================================================

localrules: all_cnvfilter

rule all_cnvfilter:
    input:
        expand("res/cnvfilter/{sample}.bed", sample = config['global']['sample-names'])
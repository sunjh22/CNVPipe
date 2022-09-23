# =================================================================================================
#     Cleaning
# =================================================================================================

def unpack_fastp_files(wildcards):
    return list(get_fastq(wildcards).values())

# output fastq must be 'trimmed', which is defined in the wrapper
rule clean_reads_se:
    input:
        sample=unpack_fastp_files
    output:
        trimmed=(
            "cleaned/{sample}.fq.gz"
            if config['settings']['keep-intermediate']['fastp']
            else temp("cleaned/{sample}.fq.gz")
        ),
        html="cleaned/{sample}-se-fastp.html",
        json="cleaned/{sample}-se-fastp.json",
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.bench.log"
    params:
        extra=config["params"]["fastp"]["se"]
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp"

rule clean_reads_pe:
    input:
        sample=unpack_fastp_files
    output:
        trimmed=(
            ["cleaned/{sample}_1.fq.gz", "cleaned/{sample}_2.fq.gz"]
            if config['settings']['keep-intermediate']['fastp']
            else temp(["cleaned/{sample}_1.fq.gz", "cleaned/{sample}_2.fq.gz"])
        ),
        html="cleaned/{sample}-pe-fastp.html",
        json="cleaned/{sample}-pe-fastp.json",
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.bench.log"
    params:
        extra=config["params"]["fastp"]["pe"]
    threads:
        config["params"]["fastp"]["threads"]
    wrapper:
        "0.64.0/bio/fastp"

rule all_fastp:
    input:
        cleaned_reads=(
            expand("cleaned/{sample}.fq.gz", sample=config["global"]["all-sample-names"])
            if is_single_end(config["global"]["all-sample-names"][0])
            else expand("cleaned/{sample}_{pair}.fq.gz", pair=[1, 2], sample=config["global"]["all-sample-names"])
        ),

localrules: all_fastp

def get_cleaned_reads(wildcards):
    if is_single_end(wildcards.sample):
        return ["cleaned/{sample}.fq.gz".format(sample=wildcards.sample)]
    else:
        return expand("cleaned/{sample}_{pair}.fq.gz", pair=[1, 2], sample=wildcards.sample)

def get_fastp_report(sample):
    if is_single_end(sample):
        return "cleaned/" + sample + "-se-fastp.json"
    else:
        return "cleaned/" + sample + "-pe-fastp.json"
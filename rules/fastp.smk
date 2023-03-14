# =================================================================================================
#     Clean reads by fastp
# =================================================================================================

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = config["global"]["samples"].loc[(wildcards.sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    else:
        return {"r1": fastqs.fq1}

def is_single_end( sample, **kargs ):
    """Return True if sample-unit is single end."""
    return pd.isnull(config["global"]["samples"].loc[(sample), "fq2"])

def unpack_fastq_files(wildcards):
    return list(get_fastq(wildcards).values())

# Since fastp wrapper is not stable, we decide to implement it by ourselves, thus use different
# commands for single-end and paired-end data.
# fastp is fast but its parallel efficiency reaches limit at 6 threads
rule clean_reads_se:
    input:
        unpack_fastq_files,
    output:
        trimmed = (
            "cleaned/{sample}.fq.gz"
            if config['settings']['keep-intermediate']['fastp']
            else temp("cleaned/{sample}.fq.gz")
        ),
        html = "cleaned/{sample}-se-fastp.html",
        json = "cleaned/{sample}-se-fastp.json",
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.bench"
    params:
        extra = config["params"]["fastp"]["se"]
    threads:
        config["params"]["fastp"]["threads"]
    conda:
        "../envs/fastp.yaml"
    shell:
        "(fastp --thread {threads} {params.extra} --in1 {input} "
        "--out1 {output.trimmed} --html {output.html} --json {output.json}) > {log} 2>&1"

rule clean_reads_pe:
    input:
        unpack_fastq_files,
    output:
        trimmed = (
            ["cleaned/{sample}_1.fq.gz", "cleaned/{sample}_2.fq.gz"]
            if config['settings']['keep-intermediate']['fastp']
            else temp(["cleaned/{sample}_1.fq.gz", "cleaned/{sample}_2.fq.gz"])
        ),
        html = "cleaned/{sample}-pe-fastp.html",
        json = "cleaned/{sample}-pe-fastp.json",
    log:
        "logs/fastp/{sample}.log"
    benchmark:
        "benchmarks/fastp/{sample}.bench"
    params:
        extra = config["params"]["fastp"]["pe"]
    threads:
        config["params"]["fastp"]["threads"]
    conda:
        "../envs/fastp.yaml"
    shell:
        "(fastp --thread {threads} {params.extra} --in1 {input[0]} --in2 {input[1]} "
        "--out1 {output.trimmed[0]} --out2 {output.trimmed[1]} "
        "--html {output.html} --json {output.json}) > {log} 2>&1"

rule multiqc_report:
    input:
        json = (
            expand("cleaned/{sample}-se-fastp.json", sample=config["global"]["sample-names"]) 
            if is_single_end(config["global"]["all-sample-names"][0])
            else expand("cleaned/{sample}-pe-fastp.json", sample = config["global"]["sample-names"])
            ),
    output:
        "cleaned/multiqc-report.html",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force -d cleaned/ -n multiqc-report -o cleaned/ -q"


localrules: all_fastp

rule all_fastp:
    input:
        cleaned_reads = (
            expand("cleaned/{sample}.fq.gz", sample=config["global"]["all-sample-names"])
            if is_single_end(config["global"]["all-sample-names"][0])
            else expand(
                "cleaned/{sample}_{pair}.fq.gz", pair=[1, 2], 
                sample=config["global"]["all-sample-names"])
        ),
        report = "cleaned/multiqc-report.html",

# Get cleaned reads no matter single-end or paired-end
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

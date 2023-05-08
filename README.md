# CNVPipe - an integrated and robust CNV calling pipeline

CNVPipe integrates several popular CNV calling methods and designed various score metrics to filter false positive callings, aims at reaching higher sensitivity and lower false discovery rate. CNVPipe uses Snakemake as backend to make it highly robust and reproducible between different machines.

## CNVPipe workflow

CNVPipe accept fastq files as input, then perform reads filtering by [fastp](https://github.com/OpenGene/fastp) and reads alignment by [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml). It will automatically determine an optimal CNV calling resolution if users not provide one. Then CNVPipe parallelly call CNVs by using five methods, including [CNVKit](https://github.com/etal/cnvkit), [CNVpytor](https://github.com/abyzovlab/CNVpytor), [cn.MOPS](http://bioconductor.org/packages/devel/bioc/html/cn.mops.html), [Smoove](https://github.com/brentp/smoove) and [Delly](https://github.com/dellytools/delly). After that, CNVPipe will adopt different strategies to recursively merge the CNV set from different methods based on the read depth. To obtain high-confidence CNVs, CNVPipe firstly use [duphold](https://github.com/brentp/duphold) to calculate the adjacent read depth of CNVs and assign "duphold score" for each CNV entry. Secondly, CNVPipe use [CNVFilteR](http://bioconductor.org/packages/release/bioc/html/CNVfilteR.html) to filter likely false positive CNVs based on B allele frequency of SNPs. In addition, CNVPipe calculate the overlap length between identified CNVs and genomic low-complexity regions to further refine the CNV set. To predict the pathogenicity of CNVs, CNVPipe compare the identified CNVs with CNVs frequently (>1%) appeared in healthy population, and assign a score based on the overlap fraction. CNVPipe also utilize [ClassifyCNV](https://github.com/Genotek/ClassifyCNV) to do a pathogenicity prediction. Finally CNVPipe will give a CNV list with various score metrics, predicted pathogenicity and possibly affected dosage-sensitive genes, and also give a list of figures showing read depth and BAF of high-confidence CNVs.

The following figure demonstrates the workflow of CNVPipe.

![CNVPipe workflow](/doc/logo/CNVPipe-workflow.png)

## Set up and usage

### Install Snakemake

You can follow the [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) to install it through [conda](https://conda.pydata.org/) or [mamba](https://github.com/mamba-org/mamba).

Following is how I install Snakemake and use it in an isolated environment.

    wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
	bash Mambaforge-Linux-x86_64.sh
	mamba create -c conda-forge -c bioconda -n snakemake snakemake
	mamba activate snakemake

### Download CNVPipe and resource files

Clone CNVPipe into your local directory

    git clone https://github.com/sunjh22/CNVPipe.git

There is a directory `resources` under CNVPipe, we recommend store all resouce files in it. Some resource files have been uploaded into GitHub, but some of them are very large, thus need to be downloaded from original source. Make separate directory for these files.

    mkdir -p resources/{refs,Delly,GATK}

Then download from following resources.
1. Download human reference genome from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/) into `resources/refs/` if no reference genome;
2. Download Delly mappbility file from [Delly](https://gear.embl.de/data/delly/) into `resources/Delly/`, remember to download all three `GRCh38` files;
3. Download GATK bundle files from [here](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/) into `resources/GATK`, all vcf files and their index are required.

### About calling resolution

A good resolution is critical for whole-genome level CNV calling, users could specifiy a value in `config` file. Based on the simulation study, we have some recommendations for different depth of data.
| Depth of data | Optimal resolution |
|:-------------:|:------------------:|
| 0.1x          | 369300             |
| 0.5x          | 73800              |
| 1x            | 36800              |
| 5x            | 7400               |
| 10x           | 3600               |
| 30x           | 1200               |

Alternatively, CNVPipe could estimate a resolution based on the median depth of input samples by ensuring each bin has 37500 base covered. To do this, we need to first run with

    snakemake --use-conda --conda-frontend mamba --conda-prefix /a/directory/where/you/store/conda/env/of/softwares --cores 10 --directory /your/working/directory autobinBydepth

This will produce a `logs/autobin/binsize.txt` file, which will be read in real CNV calling process. Actually this step will do reads filtering and alignment, thus in later CNV calling, these pre-processing steps will be skipped.

### Run CNVPipe

Notice that this workflow should always be ran under `CNVPipe` directory, users could change the working directory for different dataset. Under working directory, there should exist a folder containing fastq files of all samples, a sample sheet cantaining the location of samples' fastq files and a config file. 

1. The sample sheet could be generated by a script in `scripts/generateTable.py`, the output is `samples.tsv`. Note that [cn.MOPS](http://bioconductor.org/packages/devel/bioc/html/cn.mops.html) requires at least 6 samples to run as it applied a mixed Poission model, and [CNVKit](https://github.com/etal/cnvkit) generally requires at least one control sample. To let CNVPipe better differentiate experiment and control samples, you need to add a 'control-' prefix to the fastq file of control samples. For example, modify `myControl_1.fq.gz` to `control-myControl_1.fq.gz`. Then run the following script to generate the sample sheet.

        cd /your/working/directory
        python generateTable.py /working/directory/of/fastq/files

2. The config file could be copied from `CNVPipe` directory into working directory, users need to modify the 'absPath' and the path of various resource files in it, and specify a resolution if you do not want CNVPipe to calculate one.

3. Make sure `snakemake` environment is activated, then run CNVPipe with

        snakemake --use-conda --conda-frontend mamba --conda-prefix /a/directory/where/you/store/conda/env/of/softwares --cores 10 --directory /your/working/directory

You can use Snakemake `dry run` to do a test first

    snakemake --use-conda --conda-frontend mamba --conda-prefix /a/directory/where/you/store/conda/env/of/softwares --cores 10 --directory /your/working/directory -n

You can also run specific rules separately, for example, if you are interested in Delly, you can run with

    snakemake --use-conda --conda-frontend mamba --conda-prefix /a/directory/where/you/store/conda/env/of/softwares --cores 10 --directory /your/working/directory all_delly

Other options include:
- all_fastp: only do reads filtering
- all_bwa: only do bwa mapping and base recalibration
- all_gatk: only do GATK SNP calling
- all_freebayes: only do Freebayes SNP calling
- all_cnvkit: only run CNVKit
- all_cnvpytor: only run CNVpytor
- all_cnmops: only run cn.MOPS
- all_smoove: only run Smoove

Snakemake also provides an option to generate a graph showing all connected rules by

    snakemake --directory /your/working/directory --rulegraph | dot -Tpdf > dag.pdf

### Output of CNVPipe

The output of CNVPipe includes a summary report of reads quality, refined CNV list for each sample, and figures showing the read depth and BAF for high-confidence CNVs for each sample.

Summary report of reads quality is `cleaned/multiqc-report.html`.

CNV list is under `res/` of working directory. Results from different tools including CNVPipe, CNVKit, CNVpytor, cn.MOPS, Delly and Smoove are all under this directory. For CNVPipe, a bed file was generated for each sample, in which locations, exact copy number, various scores and pathogenicity prediction results were given. Besides, a pdf file was generated to visualize the high-confident CNVs.

## Run on clusters

To run CNVPipe on cluster environment, we need to set up some "profiles" to correctly get the resources required in our pipeline. Snakemake provided some profile templates for different cluster environments at [here](https://github.com/Snakemake-Profiles).

### Run on pbs-torque

A template profile for pbs-torque environment is provided under `profiles/pbs-torque/`. You might need to adjust the parameters such as cores and jobs in profile to fit into your cluster configuration. To run CNVPipe on pbs-torque cluster

    snakemake --use-conda --conda-frontend mamba --conda-prefix /a/directory/where/you/store/conda/env/of/softwares --profile profiles/pbs-torque --directory /your/working/directory 

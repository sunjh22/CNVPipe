# Develop a new CNV calling method - CNVPipe

Objectives of CNVPipe

1. easy to use, high sensitivity and low FDR CNV calling
2. integrate comprehensive plotting tools for visualization
3. integrate CNV annotation tools: pathogenicity prediction

## Run CNVPipe snakemake

    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data/biosoft/conda-env-cnvpipe --cores 10 --directory analysis2/

CNVPipe-token: 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

    git clone https://github.com/sunjh22/CNVPipe.git

Use bio1 server to develop with small data

Use bio2 server to run big data

## 01. Prepare test data - simulation

We want to simulate three patient samples with CNV in them and three control samples,
all downsampled to 0.01X for fast testing of our pipeline.

### Use SCNVSim to simulate SV in genome

    curl -OL https://sourceforge.net/projects/scnvsim/files/latest/download
    unzip download
    cd data/
    time java -jar ~/data/biosoft/scnvsim_1.3.1/normgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o ./simulation/simuGenome -s 0.00000001 -r 0.000000005 2>2.log &
    time java -jar ~/data/biosoft/scnvsim_1.3.1/tumorgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -k ../refs/repeatmask.txt -i ./simulation/simuGenome/normal_snvindelsim.vcf -o simulation/simuGenome -s 0.001 -l 0.001 -c 100 -a tumor1 1>1.log 2>2.log &
    $ SCNVSim Version: 1.3.1

### Use ART to simulate reads

    ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HS25 -f 1 -i simulation/simuGenome/normal_snvindelsim.fasta -l 150 -m 500 -na -o simulation/1X/control_ -p -s 10 > 1.log 2> 2.log &
    ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HS25 -f 1 -i simulation/simuGenome/tumor1_svcnvsim_clone_0.fasta -l 150 -m 500 -na -o simulation/1X/sample2_ -p -s 10 > 1.log 2> 2.log &
    gzip --fast sample_1.fq
    ART Version: v2.5.8

### Use seqkit to downsample

    cd data/simulation/1X
    parallel seqkit sample -p 0.01 -s 100 control{}_1.fq.gz -o ../../../analysis/samples-control/control{}_1.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 control{}_2.fq.gz -o ../../../analysis/samples-control/control{}_2.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 sample{}_1.fq.gz -o ../../../analysis/samples/sample{}_1.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 sample{}_2.fq.gz -o ../../../analysis/samples/sample{}_2.fq.gz ::: 1 2 3 &
    Seqkit Version: 0.10.1

`seqkit stat analysis/samples-control/*`, 116619 reads in control samples, 201149 in test samples.

## 02. Prepare snakemake config file

data: samples, control-samples, genome, refFlat

## 03. Prepare main Snakefile

Include `.smk` file step by step, testing their availability during development.

## 04. Write common.smk

Import config file, get fastq files (depending on se or pe), get patient sample and control sample names and set them to be global parameter. 
wildcard_constraints (sample names and control-sample names), valid filename and filepath.

## 05. Write pre-processing.smk

This script includes another two scripts:
- `fastp.smk` for cleaning reads, used snakemake wrapper;
- `bwamem.smk` for mapping reads using bwa-mem, also used snakemake wrapper.

## 06. Write smk for Read-depth based CNV-calling methods

### 06.1. CNVKit

Include four rules:
- `autobin` for estimating optimal bin resolution;
- `batch` for calling CNV for WGS data in batch mode;
- `segmetrics` for segmenting CNVs based on confidence interval;
- `call` for trasforming log2 depth ratio to integer copy numbers.

And a simple shell command using `awk` to extract informative columns
to another file. Columns include genomeic coordinates, integer CN,
log2 ratio, depth, probe and weight.

The refFlat and access file should be provided.

#### Principle of CNVKit

The first step is to get accessible regions. Find low mappability track: 1. get from
rqfu, in `~/data/refs/hg38/low-mappability-track/hg38.badRegions.bed`,
which includes all centromere, telomere and heterochromatin regions; 2. download from
10x genomics - black list of SV <http://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed> to
`~/data/refs/hg38/low-mappability-track/sv_blacklist.bed`, re-format it using
`grep -E -w "chr([0-9]+|[XY])" low-mappability-track/sv_blacklist.bed | sort -Vk 1 -k 2,3n > low-mappability-track/hg38.sv_blacklist.bed`.
Then use `cnvkit.py access hg38.analysisSet.fa -x hg38.badRegions.bed -x hg38.sv_blacklist.bed -o access-excludes.hg38.bed`
to get accessible regions. `cnvkit.py access` command will automatically find
sequence of N's in genome and filter regions with long N's. Then we get target bin file
by `cnvkit.py target refs/access-excludes.hg38.bed --annotate ~/data/refs/hg38/bundle/CNVKit/refFlat.txt --avg-size 20000 --split -o refs/target.bed`,
which will then be used to calculate coverage of samples in these bins.

### 06.2. cnvpytor

Include one rule with four commands:
- `rd` for extracting reads from sample bam file;
- `his` for calculating reads depth at certain resolution;
- `partition` for segmenting the bins;
- `call` for calling CNV.

The results include "CNVtype, CNVregion, CNVsize, CNVlevel, eval1, eval2, eval3, eval4, q0, pN, dG"
- `CNVlevel` - normalized read depth value, estimated integer CN value should be `round(2*CNVlevel)`;
- `eval1` - e-value (p-value multiplied by genome size divided by bin size) calculated using t-test statistics between RD statistics in the region and global - set to [0 - 0.00001];
- `eval2` - e-value from the probability of RD values within the region to be in the tails of a gaussian distribution of binned RD - set to [0 - 0.00001];
- `q0` - fraction of reads mapped with q0 quality in call region - not useful for filtering since most reads have good quality;
- `pN` - fraction of reference genome gaps (Ns) in call region - set to [0 - 0.5];
- `dG` - distance from closest large (>100bp) gap in reference genome - set to [>100kb];

Run cnvpytor

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd sample1.bam &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -his 20000 &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -partition 20000 &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -call 20000 > ../read-depth/cnvnator/sample1.2k.tsv &
    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -plot rd 20000 -o ../read-depth/cnvnator/sample1.png &

### 06.3. Control-FREEC

After reading developer's manual and testing for several times, I found it hard to include normal sample in this software,
which is really weired. So now I just need to focus on step 2 and 3. Build a config file first and run the software for
each sample.

Run freec

    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam
    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam -control analysis/mapped/control1.bam

`assess_significance.R` is not working - may need some debug.

`makeGraph.R` scirpt plot the dot plot of copy number in each window in each chromosome. gain is red, loss is blue.

### 06.4. cn.MOPS

Run cn.MOPS

    Rscript scripts/cnmops_wgs.R analysis/res/mops 20000 10 analysis/mapped/sample*.bam 2> analysis/logs/mops/call.log

Some unexpected error happened for test data, at least 6 samples are required for running mops.
I suspect it may due to the extremely low coverage of test data. The next step should be
increasing the sample number and coverage and test again. But now I have to stop at here because
Dr. Jin ask me to focus more on scCNV project. I will catch up again if I have time.
(问心无愧即可). Catch up here at Wed Oct 12 15:47:48 CST 2022.

cn.mops is not applicable to very-low-coverage data, since Windows should contain 50-100 reads each.
Three steps in cn.mops:
1. calculate read counts in windows (bins) from bam file
2. call CNV
3. calculate integer copy numbers

conda config --set channel_priority strict

## 07. Merge the results from four CNV calling tools

Emergent.

## 08. Call SNPs from data

The reason we do not use GATK for this task is that it
requires many known-variant files, which increase the
preparation load if for clinical usages. I will consider
add that option in future development.

### 08.1 freebayes

Calling sample by sample (the problem is the accuracy).
Then filter based on quality score, which is hard filter. Hope
this will not affect our purpose to correct CNV. The allele frequency
is critical here.

    freebayes -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa mapped/sample1.bam | grep 'TYPE=snp' > snps/freebayes/sample1.raw.snp.vcf

### 08.2 samtools + bcftools

Almost the same as freebayes.

## 09. Use SNPs to filter merged CNV set

### 09.1 CNVfilteR

How to install and use it in snakemake? Or I need to write a
simplified version in Python by myself. Or I need to find another
tool that was implemented in Python.

Run CNVfilteR:
1. load cnv file
2. load snv file
3. filterCNVs
4. plot single cnvs or all cnvs

CNVfilteR applies a scoring model, which is based on fuzzy logic,
to filter false positive CNVs. Briefly, copy number deletion could
be filtered if there is over 30% (default) heterozygous snvs in CNV
region; For copy number duplication, scoring model is used: the variants
with allele frequency close to 50% will be given positive value ([0,1]),
and variants with allele frequency close to 33.3% or 66.6% will be given
negative value ([-1,0]), a total score is calculated for CNV region by
summing up all the score of snvs in that CNV region; finally, positive
value indicates FP CNV and should be filtered.

This method requires as accurate SNP calling as possible. Is 10X data
really suitable for this tool? We can rank the CNVs by score but not
remove them.

.libPaths(c("/home/jhsun/R/x86_64-pc-linux-gnu-library/4.2/", "/data/jinwf/caow/R_lib/" , "/home/wangxf/R/x86_64-redhat-linux-gnu-library/3.6/"))

### 09.2 Adjacent regions

Remove all CNVs in bad genomic regions

## 10. Call CNV based on read-pair and split-read

### 10.1 Lumpy

### 10.2 Delly

## 11. Merge results from Lumpy and Delly

## 12. Combine results from read-depth based and read-pair based methods

## 13. CNV pathogenicity prediction

## 14. visualization


TODO:
1. merge some conda envs yaml files, put common tools together, only those require specific environment get a single env























## Principles of Varbin

We firstly downloaded the `k100.umap.bed` file from <https://bismap.hoffmanlab.org/>,
(Umap, Human, hg38, Single-read). This is a file showing goodzones of the genome, which
means the position can just be uniquely mapped. We then count the mappable positions in
each chromosome using `bin/varbin/mappable.py` and get `refs/k100.mappable.bed`. Then we construct
a target bin file using `bin/varbin/bin.boundary.py`, then sort this file based on chromosome using
`sort -Vk 1 -k 2,3n refs/bin.5k.test.txt > refs/bin.5k.boundary.txt`. Next, calculate the
gc content of target bins by using `bin/varbin/gc.content.py` and get `refs/bin.5k.boundary.gc.txt`.
Then we calculate the coverage of samples in target bins and correct the log2 ratio using
`bin/varbin/coverage.py`, in which `pysam` `bedcov` utils was used to calclulate the basecount in bin
regions. For example,

    time python bin/varbin/coverage.py refs/bin.5k.boundary.gc.txt data/simulation/1X/sample1.mkdup.bam analysis/read-depth/coverage/sample1.coverage.bed 1>1.log 2>2.log &


Useful commandlines:

1. sum up length of cnv region `awk 'BEGIN {s=0} {a=$3-$2;s=s+a} END {print s}' cnv.bed`

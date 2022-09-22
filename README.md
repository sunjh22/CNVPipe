# Develop a new CNV calling method - CNVPipe

Objectives of CNVPipe

1. easy to use, high sensitivity and low FDR CNV calling
2. integrate comprehensive plotting tools for visualization
3. integrate CNV annotation tools: pathogenicity prediction

## 1. Prepare test data - simulation

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

## 2. Prepare snakemake config file

data: samples, control-samples, genome, refFlat
settings:
params:

## 3. Write common.smk

import config file, get fastq files (depending on se or pe), wildcard_constraints (sample names and control-sample names), valid filename and filepath

##

Map reads, sort, mark duplicates

    time bwa mem -t 8 -M -R '@RG\tID:test0\tLB:test0\tPL:ILLUMINA\tSM:test0' ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa sample1_1.fq.gz sample1_2.fq.gz 2>2.log | samtools sort -@ 4 -o sample1.bam - 1>1.log &
    gatk MarkDuplicates -I sample1.bam -O sample1.mkdup.bam -M sample1.mkdup.metric > 1.log 2> 2.log &

Simulate 10X pair-end 100 bp data using dwgsim, map using `bowtie2`,
add read group using gatk `AddOrReplaceReadGroups`, mark duplicates
using gatk `MarkDuplicatesSpark`

Downsample 10X data to 1X using `seqkit sample -p 0.1 -s 123`, map
with `bwa mem -M`, and mark duplicates using gatk `MarkDuplicatesSpark`

    gatk AddOrReplaceReadGroups -I test0.bam -O test0.addRG.bam -LB test0 -PL ILLUMINA -PU test0 -SM test0
    gatk MarkDuplicatesSpark -I test0.addRG.bam -O test0.markdup.bam
    bowtie2 --mm -p 10 -x ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa -1 test0_1.fq.gz -2 test0_2.fq.gz | samtools sort -@ 10 > test0.bam
    seqkit sample -p 0.1 -s 123 test0_1.fq.gz -o downsample_test0_1.fq.gz
    seqkit sample -p 0.1 -s 123 test0_2.fq.gz -o downsample_test0_2.fq.gz

## Read-depth based methods

## Install `CNVKit`

 pip install numpy==1.20.0 scipy pandas matplotlib reportlab biopython pyfaidx pysam pyvcf
 pip install cnvkit

Determine read depth and bin size

 cnvkit.py autobin analysis/bam/*.bam -m wgs -b 20000 -g ~/data/refs/hg38/bundle/CNVKit/access-excludes.hg38.analysisSet.bed --annotate ~/data/refs/hg38/bundle/CNVKit/refFlat.txt

## Install `cn.MOPS`

 R
 BiocManager::install('cn.mops')

## Install `cnvnator` through conda

 conda create -n cnvnator python=3.6
 conda activate cnvnator
 conda install -c conda-forge root_base=6.20
 conda install -c conda-forge -c bioconda cnvnator

Use `cnvnator`
In the step of calculating statistics, cnvnator does not work,
throw an error "Can't find directory 'bin_1000' in file". The
problem may happen in the step of generating histogram. I cannot
resolve the problem, so I changed to another tool which used
same principle and write by the same person with Python, called
`cnvpytor.`

## Install `cnvpytor`

 git clone <https://github.com/abyzovlab/CNVpytor.git>
 cd CNVpytor
 python install --user .

Extract reads mapping

 time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd sample1.bam &

Generate histogram, 1k resolution

 time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -his 20000 &

Partition

 time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -partition 20000 &

Call CNVs

 time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -call 20000 > ../read-depth/cnvnator/sample1.2k.tsv &

Plotting

 time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -plot rd 20000 -o ../read-depth/cnvnator/sample1.png &

## Install `Control-FREEC`

 git clone <https://github.com/BoevaLab/FREEC.git>
 cd FREEC/src
 make clean
 make

Run freec

 ~/data/biosoft/FREEC/src/freec -conf analysis/read-depth/freec/config_WGS.txt -sample analysis/bam/sample1.bam -control analysis/bam-control/control1.bam

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

## Principles of CNVKit

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

Useful commandlines:

1. sum up length of cnv region `awk 'BEGIN {s=0} {a=$3-$2;s=s+a} END {print s}' cnv.bed`

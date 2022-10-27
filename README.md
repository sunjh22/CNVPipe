# Develop a new CNV calling method - CNVPipe

Objectives of CNVPipe

1. easy to use, high sensitivity and low FDR for CNV calling
2. integrate comprehensive plotting tools for visualization
3. integrate CNV annotation tools: pathogenicity prediction

## Run CNVPipe snakemake

    snakemake -np --directory analysis/ all_cnvfilter #--rerun-triggers mtime
    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data/biosoft/conda-env-cnvpipe --cores 100 --directory analysis/ all_cnvfilter
    snakemake --use-conda --conda-frontend mamba --conda-prefix /home/jhsun/data/biosoft/conda-env-cnvpipe --cores 50 --directory ~/data3/project/CNVPipe/analysis/ all_cnvfilter
    snakemake --directory analysis/ all_cnvfilter --rulegraph | dot -Tpdf > dag.pdf

CNVPipe-token: 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

    git clone https://github.com/sunjh22/CNVPipe.git

## 01. Prepare test data - simulation

For developing pipeline, we used 1X data with 6 samples and 6 controls.
For evaluating pipeline, we used 10X data with 6 samples and 6 controls.

### Use SCNVSim to simulate SV in genome

`SCNVSim` is a tool to simulate somatic SVs and CNVs in tumors, here we use it to simulate CNVs.
It has two steps: 1. simulate germline SNV and INDELs to get a normal genome, which requires 
reference genome and its length file; `-s` assigns SNV rate, `-r` assigns INDEL rate; 2. simulate 
somatic SVs and CNVs based on the first step genome and an additional repeatmask file; `-s` assigns 
structure variation rate, `-l` assigns copy neutral LOH rate, `-c` assigns number of SVs simulated, 
`-g` assigns mean copy number segment size (default 1M), `-m` assigns minimal deletion size, `-x` 
assigns maximum size of tandem duplication, `-a` prefix.

In this simulation, we need germline SNPs for later SNP calling and BAF correction, we also set
SV rate to be very low to make sure only CNV happens in the genome, so we can fully test the
performance of CNV calling tools.

`normal` genome has SNP rate 0.001, while `normal1` genome has SNP rate 0.01.
`tumor1-6` corresponds to `normal`, while `tumor7-12` corresponds to `normal1`.

    curl -OL https://sourceforge.net/projects/scnvsim/files/latest/download
    unzip download
    cd data/
    time java -jar ~/data/biosoft/scnvsim_1.3.1/normgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o ./simulation/simuGenome -s 0.01 -a normal1 1>1.log 2>2.log &
    time java -jar ~/data/biosoft/scnvsim_1.3.1/tumorgenomsim_1.3.1.jar -v ~/data/refs/hg38/analysisSet/hg38.analysisSet.size.rmMT.txt -n ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -k ~/data/refs/hg38/bundle/SCNVSim/repeatmask.txt -i ./simulation/simuGenome/normal1_snvindelsim.vcf -o simulation/simuGenome -s 0.001 -l 0.001 -c 150 -m 500 -x 10000 -a tumor7 1>1.log 2>2.log &
    $ SCNVSim Version: 1.3.1

### Use ART to simulate reads

`ART` was used to simulate paired-end sequencing reads. `-ss` model Illumina sequencing system, e.g
HSXn - HiSeqX PCR free (150bp), HS25 - HiSeq 2500 (125bp, 150bp), HSXt - HiSeqX TruSeq (150bp); `-f`
assigns coverage, `-l` assigns length, `-m` assigns mean size of DNA fragments in library; `-na` do 
not output ALN alignment file, `-s` standard deviation of DNA fragment size.

    time ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 10 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/normal1_snvindelsim.fasta -o simulation/10X/control7_ 1>1.log 2>2.log &
    time ~/data/biosoft/art_bin_MountRainier/art_illumina -ss HSXn -f 10 -l 150 -m 500 -na -p -s 10 -i simulation/simuGenome/tumor7_svcnvsim_clone_0.fasta -o simulation/10X/sample7_ >1.log 2>2.log &
    gzip -q -1 simulation/10X/control1_* &
    ART Version: v2.5.8

### Use seqkit to downsample

    cd data/simulation/1X
    parallel seqkit sample -p 0.01 -s 100 control{}_1.fq.gz -o ../../../analysis/samples-control/control{}_1.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 control{}_2.fq.gz -o ../../../analysis/samples-control/control{}_2.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 sample{}_1.fq.gz -o ../../../analysis/samples/sample{}_1.fq.gz ::: 1 2 3 &
    parallel seqkit sample -p 0.01 -s 100 sample{}_2.fq.gz -o ../../../analysis/samples/sample{}_2.fq.gz ::: 1 2 3 &
    Seqkit Version: 0.10.1

`seqkit stat analysis/samples-control/*`, 116619 reads in control samples, 201149 in test samples.

## 02. Prepare snakemake main file and config file

Main file: include `.smk` file step by step, testing their availability during development.
Config file: data: samples, control-samples, genome, refFlat

## 03. Write common.smk

Import config file, get fastq files (depending on se or pe), get patient sample and 
control sample names and set them to be global parameter. wildcard_constraints (sample 
names and control-sample names), valid filename and filepath.

## 04. Write pre-processing.smk

This script includes another two scripts:
- `fastp.smk` for cleaning reads, use snakemake wrapper;
- `bwamem.smk` for mapping reads using bwa-mem, also use snakemake wrapper.

## 05. Write smk for Read-depth based CNV-calling methods

We need to be clear about some features of these methods:
1. require normal samples? if yes, how many?
2. call germline CNV sample by sample or in batch
3. exclude low-mappability region? require reference file?

### 05.1. CNVKit

Require normal samples (no limit on the number), reference genome and exclude low-map region.
Call germline CNV in batch.

Include four rules:
- `autobin` for estimating optimal bin resolution;
- `batch` for calling CNV for WGS data in batch mode;
- `segmetrics` for segmenting CNVs based on confidence interval;
- `call` for transforming log2 depth ratio to integer copy numbers.

And a simple shell command using `awk` to extract informative columns
into another file. Columns include genomeic coordinates, integer CN,
log2 ratio, depth, probe and weight.

The refFlat and access file should be provided for annotation and
discarding reads in low-mappbility region.

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

### 05.2. cnvpytor

Do not require normal sample and reference genome and do not exclude low-map region.
This method solely depends on comparing read depth between adjacent region (?) and regions
with similar GC content. Call germline CNV sample by sample.

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

Get read depth signal

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -chrom $(seq -f 'chr%g' 1 22) chrX chrY -rd sample1.bam &

Assign bin size, calculate read depth histogram

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -his 20000 &

Segment by mean-shift method

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -partition 20000 &

Call CNV

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -call 20000 > ../read-depth/cnvnator/sample1.2k.tsv &

Plot

    time cnvpytor -root ../read-depth/cnvnator/sample1.pytor -plot rd 20000 -o ../read-depth/cnvnator/sample1.png &

Finally use a python script to extract evalue, pN and dG into `info` column

`cnvpytor` accepts SNP data and call CNV based on both RD and BAF, but we want to separate them,
just use its RD function to call CNV, and filter by BAF after merging CNV from other tools.

#### Principle of CNVpytor

Firstly, extract RD signal in 100-bp intervals (pysam), get bin-level RD signal, perform GC biase
correction (100-bp GC content is pre-calculated), then merge and partition by **mean-shift**, CNV
call follows, all requires statistical test of RD between CNV segments and whole-genomic regions.

Filter out HET SNPs that are inaccessible to short read, defined by strict mask from 1000G project,
BAF likelihood of bin is calculated by multiplying the likelihood of single SNPs in the bin (which
is modelled by a beta distribution ?).

Annotating CNV with genes (inside, covering or intersecting left/right breakpoints)

Merge calls over multiple regions: >50% reciprocal overlap

### 05.3. Control-FREEC

`Control-FREEC` is written in C++.

This tool can be run in two modes: 1. for WGS data without normal reference, a directory with fasta
of single chromosomes or a GC-content profile need to be provided; 2. for WES data, a normal 
reference must be provided. For mode1, a `GEM` mappability file could be provided to remove bad
genomic regions. In my testing, control parameter is not working even if you provide it, do not know 
the reason.

Generate GC content profile for specific size of window. In this case, `chrFiles` should be set to
the directory with fasta sequences of single chromosomes, and `window` set to specific size.

    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam

Run freec with specified GC-content profile and mappability file

    ~/data/biosoft/FREEC/src/freec -conf analysis/temp/freec/config_WGS.txt -sample analysis/mapped/sample1.bam

After checking the results, we found that even a mappability file is provided, some telomere regions 
are called to be CNV, which either results from the problem of mappability file or from the tool
itself. At least, this result makes its all calling result not reliable.

#### Principle of Control-FREEC

Control-FREEC takes as an input aligned reads, then constructs and normalizes the copy number 
profile (by GC profile or normal read coverage), constructs the B-allele frequency (BAF) profile, 
segments both profiles (Lasso-based algorithm), ascribes the genotype status to each segment using 
both copy number and allelic frequency information, then annotates genomic alterations.

### 05.4. cn.MOPS

Require at least 6 normal samples, do not require reference genome and do not exclude low-map 
region. Call CNV in batch. This method call CNV by comparing read depth between patient sample and 
control sample.

cn.mops is not applicable to very-low-coverage data, since Windows should contain 50-100 reads each.
Three steps in cn.mops:
1. calculate read counts in windows (bins) from bam file
2. call CNV
3. calculate integer copy numbers

Run cn.MOPS

    Rscript scripts/cnmops_wgs.R analysis/res/mops 20000 10 analysis/mapped/sample*.bam 2> analysis/logs/mops/call.log

Some unexpected error happened for test data, at least 6 samples are required for running mops.
I suspect it may due to the extremely low coverage of test data. The next step should be
increasing the sample number and coverage and test again. But now I have to stop at here because
Dr. Jin ask me to focus more on scCNV project. I will catch up again if I have time.
(问心无愧即可). Catch up here at Wed Oct 12 15:47:48 CST 2022.

## 06. Call CNV based on read-pair and split-read

### 06.1 Lumpy

Lumpy is written in C.
Five steps to call CNV by Lumpy. 1. get discordant reads (insert size over expected size or
mate pair mapped to different chromosomes); 2. get split reads by Lumpy script; 3. call CNV by
lumpyexpress; 4. genotype CNV by svtypes, a quality score was given (in QUAL filed of vcf file),
but how to select CNV based on score need further consideration; 5. extract CNV.

    samtools view -b -F 1294 analysis/mapped/sample1.bam | samtools sort -@ 8 -o analysis/temp/lumpy/sample1.discordant.bam -
    samtools view -h analysis/mapped/sample1.bam | extractSplitReads_BwaMem -i stdin | samtools sort -@ 8 -o analysis/temp/lumpy/sample1.split.bam -
    lumpyexpress -B analysis/mapped/sample1.bam -S analysis/temp/lumpy/sample1.split.bam -D analysis/temp/lumpy/sample1.discordant.bam -o analysis/temp/lumpy/sample1.vcf
    svtyper -i analysis/temp/lumpy/sample1.vcf -B analysis/mapped/sample1.bam -l analysis/temp/lumpy/sample1.bam.json > analysis/temp/lumpy/sample1.gt.vcf
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL\n' analysis/temp/lumpy/sample1.gt.vcf | egrep 'DUP|DEL' > analysis/res/lumpy/sample1.bed

A coverage calculation script was provided `scripts/get_coverages.py` for later reference.

A low mappbility track file provided by Heng Li is suggested, I downloaded and extracted
the file at Mon Oct 17 08:59:27 +08 2022, the file is in 
`~/refs/hg38/low-mappability-track/HengLi-lowMapTrach/hg38.exclude.hengli.bed`.
Overall length of low-map region is 228,061,378.

### 06.2 Smoove

`smoove` is a Lumpy wrapper, which is easier to use and faster.

An bad region file was downloaded according to the author of smoove. Overall length is 119,556,880,
which is almost half shorter than Heng Li's. But when we remove non-cannonical chromosomes, its 
length became 3,042,685, which is definitely not correct. 
The file is in `low-mappability-track/smoove/hg38.exclude.bed`.

    wget https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed

Our own bad region file is in `low-mappability-track/hg38.badRegions.bed`, which includes 
centromere, telomere and heterochromatin, the length is 199,765,358.

Another low-map track is downloaded from 10x genomics in 
`low-mappability-track/10x/hg38.sv_blacklist.bed`, its overall length is 212,765,070.
I think this one might be the best one to use.

Test `smoove` calling

    smoove call --outdir analysis/temp/smoove/ --exclude ~/data/refs/hg38/exclude.cnvnator_100bp.GRCh38.20170403.bed --name sample4 --fasta ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -p 1 --genotype analysis/mapped/sample4.bam

After comparing the results from Lumpy and Smoove, we found no big difference, so we will use smoove 
to replace lumpy in our CNVPipe.

Use `duphold` to genotype. Three additional fields will be added into `FORMAT` per sample:
DHFC - "fold-change for the variant depth relative to the rest of the chromosome the variant was found on"
DHFFC - "fold-change for the variant depth relative to Flanking regions"
DHBFC - "fold-change for the variant depth relative to bins in the genome with similar GC-content"

`duphold` could be directly downloaded as an exectuable binary file.

    smoove genotype -d -x -p 1 --name sample4 --outdir analysis/temp/smoove-genotype/ --fasta ~/data3/refs/hg38/analysisSet/hg38.analysisSet.fa --vcf analysis/temp/smoove/sample4-smoove.genotyped.vcf.gz analysis/mapped/sample4.bam

    Version: Smoove 0.2.8

The authors suggest use DHFFC<0.7 to filter deletion and DHBFC>1.3 to filter duplication

    bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3)' analysis/temp/smoove-genotype/sample4-smoove.genotyped.vcf.gz

Use `bcftools` to extract informative columns

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%QUAL[\t%DHFC\t%DHFFC\t%DHBFC]\n' analysis/temp/smoove-genotype/sample4-smoove.genotyped.vcf.gz

Follow author's suggestion to filter CNV based on DHFFC and DHBFC by `awk`

    awk '$4=="DEL" && $7<0.7 {print $0} $4=="DUP" && $8>1.3 {print $0}' analysis/res/smoove/sample1.bed | wcl

`duphold` also support SNP vcf file for filtering, just like `CNVfilteR`.

### 06.3 Delly

Delly is written in C++. It is primarily designed for calling structural variant, both germline and 
somatic (cancers). It also supports CNV calling now.

Mappbility files are mandantory for Delly to call CNV, we downloaded map files
from [here]<https://gear.embl.de/data/delly/> at Sat Oct 15 17:08:30 +08 2022. All three files are 
required.

Delly by default divide genome into 10kb-mappable bins, but we can set window size by `-i`
parameter. `-u` is for segmentation, `-l` is for using delly SV calling to refine breakpoints.

So we try to call SV first, and use its output as input for CNV calling. 
`delly call` only use paired-end and split-read information, read-depth is not used.

    delly call -g ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -x ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed -o analysis/temp/delly/sample1.sv.bcf analysis/mapped/sample1.bam

Use `delly cnv` to call CNV

    delly cnv -u -i 20000 -g ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -l analysis/temp/delly/sample1.sv.bcf -m ~/data/refs/hg38/bundle/delly/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz -c analysis/temp/delly/sample1.cov.gz -o analysis/temp/delly/sample1.cnv.bcf analysis/mapped/sample1.bam

Use `duphold` to do genotyping

    duphold -t 4 -v analysis/temp/delly/sample4.filtered.bcf -b analysis/mapped/sample4.bam -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa -o analysis/temp/delly/sample4.duphold.vcf

`delly classify` will select CNVs with 'PASS' (filter low-quality CNVs), bad thing is that it
transforms all quality score into 10000, so we decided to filter CNVs by grep and keep
the original quality score. One thing need to be noticed is that in some cases, value in `FILTER`
field is not matched with `FT` tag in `format` field. One could be PASS while the other one
could be LowQual. We just use the `FILTER`=="PASS" to select all good-quality CNVs. (strict mode)

Use `bcftools` to extract informative columns

    bcftools query -f '%FILTER\t%CHROM\t%POS\t%INFO/END[\t%RDCN\t%DHFC\t%DHFFC\t%DHBFC]\n' analysis/temp/delly/sample4.duphold.vcf | grep 'PASS' | cut -f 2- > analysis/res/delly/sample4.bed

In `INFO` field, `MP` means fraction of mappable positions; `GCF` means GC fraction (added by 
duphold); these should not be used for CNV filtering.

We also need to filter CNV based on DHFFC and DHBFC as stated in `smoove` section.

    awk '$4<2 && $7<0.7 {print $0} $4>2 && $8>1.3 {print $0}' analysis/res/delly/sample1.bed | wcl

Some R scripts were provided for plotting CNV, which could be a reference for later usage.

    Rscript ~/data/biosoft/delly/R/rd.R analysis/temp/delly/sample1.cov.gz analysis/temp/delly/sample1.segment.bed

    Version: Delly 1.1.5

## 07. Call SNPs from data

The reason we do not use GATK for this task is that it
requires many known-variant files, which increase the
preparation load if for clinical usages. I will consider
add that option in future development.

### 07.1 freebayes

Call SNP sample by sample (the problem is the accuracy). Then filter based on quality score
and read depth, this is a critical step. We should set different threshold for different depth
data. The accuracy of SNP calling will affect the correction of CNV in later step.

    freebayes -f ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa mapped/sample1.bam > snps/freebayes/sample1.raw.vcf
    bcftools filter -O v -o snps/freebayes/sample1.filtered.vcf -s LOWQUAL -e 'QUAL<10 || FMT/DP <5' --SnpGap 5 --set-GTs . snps/freebayes/sample1.raw.vcf
    bcftools view -v snps snps/freebayes/sample1.filtered.vcf > snps/freebayes/sample1.snp.vcf

### 07.2 samtools + bcftools

Almost the same as freebayes.

### 07.3 GATK

Download GATK bundle resource file from [here]<https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/>
at Thu Oct 27 11:29:27 +08 2022. I suspect bundle files in our server might have broken or I misuse
some of files in workflow. After downloading and checking with md5sum, I found they are totally
identical.

When runing `variantRecalibrator`, an error happened: "Bad input: Found annotations with zero 
variance". They must be excluded before proceeding. This might be due to the very-low-coverage of
our testing data.
Another error is "Positive training model failed to converge.  One or more annotations (usually MQ) 
may have insufficient variance", this is caused by too less number of SNPs identified in GATK 
`HaplotypeCaller`. We need to consider re-simulate sequencing data with high SNP rate.

## 08. Merge the results from all CNV calling tools

    python ~/data/project/CNVPipe/scripts/mergeCNV.py analysis/res/cnvkit/sample1.bed analysis/res/cnvpytor/sample1.bed analysis/res/freec/sample1.bed analysis/res/mops/sample1.bed analysis/res/smoove/sample1.bed ~/data/refs/hg38/low-mappability-track/10x/hg38.sv_blacklist.bed analysis/res/delly/sample1.bed sample1.merge.bed >merge.log

Except merging the results from four calling tools, we try to remove genomic bad regions 
(centromere, telomere and heterochromatin). The bad region file is in 
`~/data/refs/hg38/low-mappability-track/10x//hg38.sv_blacklist.bed`. Besides, we need to make sure the 
merged file adapt to the input format of `CNVfilteR`.

We tried to merge CNV results from all tools including cnvkit, cnvpytor, freec, mops, smoove 
and delly.

? What is heterochromatin

## 09. Score by overlap fractions with genomic bad regions

The basic idea is that the higher overlap fraction, the lower score CNV region will get.

## 10. Score by duphold results

Following last step, we can transform the bed file into vcf file as `duphold` requires VCF format
as input. Then we use it to calculate adjacent read depth for CNV region, and design a strategy
to score CNV regions. Transforming script is in `scripts/bed2vcf.py`.

    python ~/data/project/CNVPipe/scripts/bed2vcf.py sample1.bed ~/data/refs/hg38/analysisSet/hg38.analysisSet.fa.fai sample1.vcf

## 11. Score by BAF-correction

`CNVfilteR`

How to install and use it in snakemake? Or I need to write a simplified version in Python by myself. 
Or I need to find another tool that was implemented in Python. Directly download in R script if not 
installed.

Run CNVfilteR:
1. load cnv file
2. load snv file
3. filterCNVs
4. plot single cnvs or all cnvs

CNVfilteR applies a scoring model, which is based on fuzzy logic, to filter false positive CNVs. 
Briefly, copy number deletion could be filtered if there is over 30% (default) heterozygous snvs in 
CNV region; For copy number duplication, scoring model is used: the variants with allele frequency 
close to 50% will be given positive value ([0,1]), and variants with allele frequency close to 33.3% 
or 66.6% will be given negative value ([-1,0]), a total score is calculated for CNV region by 
summing up all the score of snvs in that CNV region; finally, positive value indicates FP CNV and 
should be filtered.

This method requires as accurate SNP calling as possible. Is 10X data really suitable for this tool? 
We can rank the CNVs by score but not remove them.

The script is in `scripts/cnvFilteR.R`

## 12. Integrate results from different tools


## 13. CNV pathogenicity prediction

### 13.1 ClassifyCNV

ClassifyCNV cannot be downloaded from bioconda, and its result directory is weired,
we need to figure out how to make it work in our CNVPipe.

### 13.2 X-CNV

XCNV only supports hg19, which is outdated, so we give up using this tool to predict
CNV pathogenicity.

## 14. visualization
























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
